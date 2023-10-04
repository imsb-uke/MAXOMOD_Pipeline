1. About NESmapper

NESmapper is a perl program to predict CRM1-dependent, leucine-rich nuclear export signals (NESs) in proteins or peptides. NESmapper calculates an NES score of a peptide consisting of 11~14 amino acids using scores in NES profiles. NES profiles are matrices of score for four classes of NES (class 1a/3, class 1b, class 1c, and class 2), in which each score represents a relative NES activity level corresponding to the position (column) and amino acid (line) of a template NES indicated at the top line. The scores indicated in the profiles have been determined by assaying nuclear export activities of NES mutants of the template NES in mammalian cells. Two files containing NES profiles are provided in this package, NES_profiles.txt and NES_profiles_trained.txt, the latter of which has been trained with experimentally confirmed NES data to optimize the scores.

2. Prerequisites

Perl 5.6 or laterIf your perl is not in the standard location (/usr/bin/perl), edit the first line of the Coval scripts or create a link of your perl executable to /usr/bin.

3. Usage

Set the environmental variable PATH to the directory containing NESmapper.pl and NES_profile_trained.txt or use {$script_path: absolute_path_to_NESmapper_directry}/NESmapper.pl

{$script_path}/NESmapper.pl [options] <input_fasta> or <STDIN>

The input file is a fasta file containing protein sequence(s). Unless an input file is given, you can paste a single amino acid sequence onto the command line following a message 'Paste your amino acid sequence' and should press the return key one or two times. By default, NES_profiles_trained.txt located in the $script_path is used for NES profiles and the number of predicted NESs is only output on the console. If the details of the predicted NESs are required, specify -o option with an output file name.

<Available options>

  --minscore or -ms <INT>   minimum score to judge whether the predicted NESs are true or false [default: 2]
  --disable_out or -d <STR> disable to output predicted NESs (only NESs with >= minscore are output by default) [default: false (output)]
  --profile or -p <STR>     NES profile file [default: NES_profiles_trained.txt]
  --help or -h              output help message

4. Other

An accompanied script, NESmapper_ValidNES.pl, is a program to predict NESs from the data from the ValidNES database and to evaluate the result. Please see the usage by conducting NESmapper_ValidNES.pl -h.


5.LicenseGMcloserCopyright (c) 2014, Shunichi Kosugi, Kazusa DNA Research InstituteLicensed under the GNU General Public License Version 3 license (http://www.gnu.org/licenses/lgpl.txt).6. Contact

Shunichi Kosugi  Kazusa DNA Research InstituteEmail: skosugi@kazusa.or.jp