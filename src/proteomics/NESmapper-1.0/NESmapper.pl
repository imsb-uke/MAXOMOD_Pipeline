#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

my $NES_profile = "$Bin/NES_profiles_trained.txt";
my $threshold_score = 2;
my $cut_off = -10;
my $notout_predicted_NES = 0;
my $help;

print STDERR "$Bin\n";

GetOptions(
    'profile|p=s' => \$NES_profile,
    'minscore|ms=i' => \$threshold_score,
    'disable_out|d' => \$notout_predicted_NES,
    'help|h' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  NES_mapper.pl [options] <input_fasta> or <STDIN>
  Options:
   --profile or -p      NES profile file
   --minscore or -ms    minimum score to judge whether the predicted NESs are true or false [default: 2]
   --disable_out or -d  disable to output predicted NESs [default: false]
   --help or -h         output help message

=cut

if (!-e $NES_profile){
    $NES_profile = "$Bin/$NES_profile";
}

# for the standard 14-aa NES profiles, scan 14aa-window

my $file = '';
my %input_seq;


my $max_hydrophobicity_rate = 0.4;      # allowable rate of hydrophobic residues in NES spacer region

if (@ARGV > 0){
    $file = shift (@ARGV);
}
else{
    print STDERR "Paste your amino acid sequence.\n";
    while (my $line = <>) {
        if ($line =~ /[bjouxz]/i){
            print "Error in sequence; your sequence contains an illegalcharacter(s).\n";
            exit;
        }
        elsif ($line eq "\n"){
            last;
        }
        else {
            $line =~ s/[\W*\d*\b*\r*\s*]//g;
            my $text_seq .= uc $line;
            if (length $text_seq < 14){
                my $added_num = 14 - length $text_seq;
                $text_seq = ('S' x $added_num) . $text_seq;
            }
            $input_seq{'test_nes'} = $text_seq;
        }
    }
}

my $Seq = '';
my $header = '';
my %protein_id;
my %common_protein;

if ($file ne ''){
    open (FILE, $file) or die "$file: $!";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^$/);
        if ($line =~ /^>(.+)/){
            if ($Seq ne ''){
                if (length $Seq < 14){
                    my $added_num = 14 - length $Seq;
                    $Seq = ('S' x $added_num) . $Seq;
                }
                $input_seq{$header} = $Seq;
                $Seq = '';
            }
            my @header = split (/\s+/, $1);
            $header = $header[0];
            $protein_id{$header} = $header;
            if (($header =~ /^S\d+$/) and (@header > 2)){   # for ValidNES_data_LMBtested.fasta
                $protein_id{$header} = $header[2];
            }
        }
        else{
            $line =~ s/[\W*\d*\b*\r*\s*]//g;
            $Seq .= uc $line;
        }
    }
    $input_seq{$header} = $Seq;
    $Seq = '';
}


my %profile_1a;
my %profile_1b;
my %profile_1c;
my %profile_2;
my @score_line;
my $NES_class = '';
my @aa_line = ('B', 'L', 'I', 'V', 'M', 'C', 'F', 'Y', 'W', 'T', 'P', 'A', 'G', 'S', 'H', 'R', 'K', 'Q', 'N', 'D', 'E');
open (FILE, $NES_profile) or die "$!";
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^$/);
    if ($line =~ /^(class\S+)/){
        if (@score_line > 0){
            my $count_line = 0;
            foreach my $sline (@score_line){
                $count_line ++;
                my @scores = split (/\t/, $sline);
                my $count_pos = 0;
                foreach (@scores){
                    if ($count_pos > 0){
                        ${$profile_1a{$aa_line[$count_line]}}{$count_pos} = $_ - 8 if ($NES_class eq 'class-1a');
                        ${$profile_1b{$aa_line[$count_line]}}{$count_pos} = $_ - 4 if ($NES_class eq 'class-1b');
                        ${$profile_1c{$aa_line[$count_line]}}{$count_pos} = $_ - 8 if ($NES_class eq 'class-1c');
                        ${$profile_2{$aa_line[$count_line]}}{$count_pos} = $_ - 3 if ($NES_class eq 'class-2');
                    }
                    $count_pos ++;
                }
            }
            @score_line = ();
        }
        $NES_class = $1;
    }
    else{
        push @score_line, $line if ($line =~ /^[A-Z]/);
    }
}
close (FILE);

# Calculate activity scores of CRM1-dpendent NESs consisting of 14 consecutive amino acids within the query sequence by scanning from the N-terminus to the C-terminus.

my $sum_1a = 0;
my $sum_1b = 0;
my $sum_1c = 0;
my $sum_2 = 0;
my $total_query = 0;
my %sum_all;
my %duplicate;
my %hit_query;
my %hit_nes;

if (($notout_predicted_NES == 0) or ($file eq '')){
    print "<query>\t<pos>\t<nes>\t<score>\t<nes-class>\n" if ($notout_predicted_NES == 0);
    print "<query>\t<pos>\t<nes>\t<score>\t<nes-class>\n" if ($file eq '');
}
foreach my $name (sort keys %input_seq){
#    next if (exists $duplicate{$protein_id{$name}});
    if ($name =~ /^S\d+$/){                                         # for ValidNES_data_LMBtested.fasta
        $duplicate{$protein_id{$name}} ++;
    }
    $total_query ++;
    my @sum_array;
    my %hit_pos;
    my %hit_pos_2;
    my %hit_pos_3;
    my $seq = $input_seq{$name};
    for (my $i = 0; $i < 3; $i++){      # Calculate scores of NES beginning at the position 4.
        $sum_1a = 0;
        $sum_1b = 0;
        $sum_1c = 0;
        $sum_2 = 0;
        my @spacer_residue_1a;
        my @spacer_residue_1b;
        my @spacer_residue_1c;
        my @spacer_residue_2;
        my $nes_sub = substr ($seq, 0, $i + 11);
        my @extracted_residues = split (//, $nes_sub);
        push @spacer_residue_1a, $extracted_residues[$i + 1], $extracted_residues[$i + 2], $extracted_residues[$i + 3], $extracted_residues[$i + 5], $extracted_residues[$i + 6], $extracted_residues[$i + 8];
        push @spacer_residue_1b, $extracted_residues[$i + 2], $extracted_residues[$i + 3], $extracted_residues[$i + 5], $extracted_residues[$i + 6], $extracted_residues[$i + 8];
        push @spacer_residue_1c, $extracted_residues[$i], $extracted_residues[$i + 1], $extracted_residues[$i + 2], $extracted_residues[$i + 4], $extracted_residues[$i + 5], $extracted_residues[$i + 6], $extracted_residues[$i + 8] if ($i > 0);
        push @spacer_residue_2, $extracted_residues[$i + 3], $extracted_residues[$i + 5], $extracted_residues[$i + 6], $extracted_residues[$i + 8];
        my @HH_1a = grep {/[LIMVF]/} @spacer_residue_1a;
        my @HH_1b = grep {/[LIMVF]/} @spacer_residue_1b;
        my @HH_1c = grep {/[LIMVF]/} @spacer_residue_1c if ($i > 0);
        my @HH_2 = grep {/[LIMVF]/} @spacer_residue_2;
        $sum_1a -= int (@HH_1a / 6 * 7) if (@HH_1a / 6 >= $max_hydrophobicity_rate);
        $sum_1b -= int (@HH_1b / 5 * 7) if (@HH_1b / 5 >= $max_hydrophobicity_rate);
        $sum_1c -= int (@HH_1c / 7 * 7) if (@HH_1c / 7 >= $max_hydrophobicity_rate);
        $sum_2 -= int (@HH_2 / 4 * 7) if (@HH_2 / 4 > $max_hydrophobicity_rate);
        for (my $j = 0; $j < $i + 11; $j++){
            my $pf_pos = $j - $i + 4;
            $sum_1a += ${$profile_1a{$extracted_residues[$j]}}{$pf_pos};
            $sum_1b += ${$profile_1b{$extracted_residues[$j]}}{$pf_pos};
            $sum_1c += ${$profile_1c{$extracted_residues[$j]}}{$pf_pos} if ($i > 0);
            $sum_2 += ${$profile_2{$extracted_residues[$j]}}{$pf_pos};
        }
        my $extracted_residues = join ('', @extracted_residues);
        $extracted_residues = ('S' x (3 - $i)) . $extracted_residues;
        push @sum_array, 4 - $i, $extracted_residues, $sum_1a + 8, 'class_1a', $name if ($sum_1a + 8 >= $cut_off);
        push @sum_array, 4 - $i, $extracted_residues, $sum_1b + 4, 'class_1b', $name if ($sum_1b + 4 >= $cut_off);
        push @sum_array, 4 - $i, $extracted_residues, $sum_1c + 8, 'class_1c', $name if ($sum_1c + 8 >= $cut_off) and ($i > 0);
        push @sum_array, 4 - $i, $extracted_residues, $sum_2 + 3, 'class_2', $name if ($sum_2 + 3 >= $cut_off);
    }
    
    for (my $i = 0; $i <= length ($seq) - 14; $i++){        # Calculate scores of full-length NES.
        $sum_1a = 0;
        $sum_1b = 0;
        $sum_1c = 0;
        $sum_2 = 0;
        my @spacer_residue_1a;
        my @spacer_residue_1b;
        my @spacer_residue_1c;
        my @spacer_residue_2;
        my $nes_sub = substr ($seq, $i, 14);
        my @extracted_residues = split (//, $nes_sub);
        push @spacer_residue_1a, $extracted_residues[4], $extracted_residues[5], $extracted_residues[6], $extracted_residues[8], $extracted_residues[9], $extracted_residues[11];
        push @spacer_residue_1b, $extracted_residues[5], $extracted_residues[6], $extracted_residues[8], $extracted_residues[9], $extracted_residues[11];
        push @spacer_residue_1c, $extracted_residues[3], $extracted_residues[4], $extracted_residues[5], $extracted_residues[7], $extracted_residues[8], $extracted_residues[9], $extracted_residues[11];
        push @spacer_residue_2, $extracted_residues[6], $extracted_residues[8], $extracted_residues[9], $extracted_residues[11];
        my @HH_1a = grep {/[LIMVF]/} @spacer_residue_1a;
        my @HH_1b = grep {/[LIMVF]/} @spacer_residue_1b;
        my @HH_1c = grep {/[LIMVF]/} @spacer_residue_1c;
        my @HH_2 = grep {/[LIMVF]/} @spacer_residue_2;
        $sum_1a -= int (@HH_1a / 6 * 7) if (@HH_1a / 6 >= $max_hydrophobicity_rate);
        $sum_1b -= int (@HH_1b / 5 * 7) if (@HH_1b / 5 >= $max_hydrophobicity_rate);
        $sum_1c -= int (@HH_1c / 7 * 7) if (@HH_1c / 7 >= $max_hydrophobicity_rate);
        $sum_2 -= int (@HH_2 / 4 * 7) if (@HH_2 / 4 > $max_hydrophobicity_rate);
        for (my $j = 0; $j < 14; $j++){
            my $pf_pos = $j + 1;
            $sum_1a += ${$profile_1a{$extracted_residues[$j]}}{$pf_pos};
            $sum_1b += ${$profile_1b{$extracted_residues[$j]}}{$pf_pos};
            $sum_1c += ${$profile_1c{$extracted_residues[$j]}}{$pf_pos};
            $sum_2 += ${$profile_2{$extracted_residues[$j]}}{$pf_pos};
        }
        my $extracted_residues = join ('', @extracted_residues);
        push @sum_array, $i + 1, $extracted_residues, $sum_1a + 8, 'class_1a', $name if ($sum_1a + 8 >= $cut_off);
        push @sum_array, $i + 1, $extracted_residues, $sum_1b + 4, 'class_1b', $name if ($sum_1b + 4 >= $cut_off);
        push @sum_array, $i + 1, $extracted_residues, $sum_1c + 8, 'class_1c', $name if ($sum_1c + 8 >= $cut_off);
        push @sum_array, $i + 1, $extracted_residues, $sum_2 + 3, 'class_2', $name if ($sum_2 + 3 >= $cut_off);
    }
    my $count = 0;
    while (@sum_array > 4){
        my ($pos, $nes, $score, $class, $id) = splice (@sum_array, 0, 5);
        my $name_pos = "$name-$pos";
        if (length $seq > 25){
            my $up_seq = '';
            my $dw_seq = '';
            my $flank_score = 1;
            my $hh_pol_r = 0;
            my $charge = -100;
            if ($pos - 26 >= 0){
                $up_seq = substr ($seq, $pos - 26, 25);
            }
            elsif ($pos >= 5){
                $up_seq = substr ($seq, 0, $pos - 1);
            }
            if ($pos + length ($nes) - 1 <= length ($seq) - 25){
                $dw_seq = substr ($seq, $pos + length ($nes) - 2, 25);
            }
            elsif ($pos + length ($nes) - 1 <= length ($seq) - 4){
                $dw_seq = substr ($seq, $pos + length ($nes) - 2);
            }
            if ($up_seq ne ''){
                my $hh = 0;
                my $pol = 0;
                $hh = $up_seq =~ s/[LMVIFCW]/X/g;
                $pol = $up_seq =~ s/[SKRHDENQP]/X/g;
                $pol = 1 if ($pol eq '') or (!defined $pol);
                $hh_pol_r = int ($hh / $pol * 1000) / 10;
            }
            if ($dw_seq ne ''){
                my $acid = 0;
                my $base = 0;
                $acid = $dw_seq =~ s/[DE]/X/g;
                $base = $dw_seq =~ s/[KRH]/X/g;
                $charge -= $acid;
                $charge += $base;
            }
            if (($hh_pol_r > 0) and ($hh_pol_r <= 30)){
                $flank_score *= 2.5;
            }
            elsif ($hh_pol_r <= 40){
                $flank_score *= 2;
            }
            elsif ($hh_pol_r <= 50){
                $flank_score *= 1.4;
            }
            elsif ($hh_pol_r <= 60){
                $flank_score *= 1;
            }
            elsif ($hh_pol_r <= 80){
                $flank_score *= 0.6;
            }
            elsif ($hh_pol_r > 80){
                $flank_score *= 0.5;
            }
            if ($charge <= -4){
                $flank_score *= 1.8;
            }
            elsif ($charge <= 0){
                $flank_score *= 1;
            }
            else{
                $flank_score *= 0.6;
            }
            $score *= $flank_score if ($score > 0);
        }
        if ($score >= $threshold_score){
            next if (exists $hit_pos{$pos});
            $hit_query{$name} = 1;
            $hit_nes{$name_pos} = 1;
            for (my $i = 0; $i < 14; $i++){
                my $pos_2 = $pos + $i;
                $hit_pos{$pos_2} = 1;
            }
            my $pos_05d = sprintf ("%05d", $pos);
            my $id2 = "$id-$pos_05d-$nes";
            if (($notout_predicted_NES == 0) or ($file eq '')){
                print "$name\t$pos\t$nes\t$score\t$class\n";
            }
        }
    }
}
if ($file eq ''){
    print "\n";
}

my $hit_query_num = scalar keys %hit_query;
my $hit_nes_num = scalar keys %hit_nes;
my $rate_hit = int ($hit_query_num / $total_query * 1000) / 10;

print "Query number: $total_query\n";
print "Hit query number (score >= $threshold_score): $hit_query_num ($rate_hit%)\n";
print "Hit NES number (score >= $threshold_score): $hit_nes_num\n";
