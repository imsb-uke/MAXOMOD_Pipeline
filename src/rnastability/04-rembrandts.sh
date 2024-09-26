#!/bin/bash
set -e

####################### define executables
deseq="./src/rnastability/rembrandts/DESeq.R"
rembrandts="./src/rnastability/rembrandts/REMBRANDTS.R"

####################### identify the input arguments
jobid=$$
metadata=$1
refdir=$2
stringency=$3
fitmode=$4
outdir=$5

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash REMBRANDTS.sh <metadata.txt> <inputDir> <stringency> <biasMode> <outDir>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input metadata file: "$metadata
echo "Reference directory for HTSeq-Count files: "$refdir
echo "Stringency for filtering measurements: "$stringency
echo "The mode of identifying bias parameters: "$fitmode
echo "Ouput Directory: "$outdir

if [ -e "$metadata" ]; then
	echo "Metadata file found."
else
	echo "ERROR: Metadata file was not found."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder

####################### run DESeq
Rscript $deseq $jobid $metadata $refdir


####################### define output path
out_folder="./out/"$jobid
mkdir -p $out_folder
mkdir -p $out_folder"/sampleScatterplots"

####################### run REMBRANDTS

Rscript $rembrandts $jobid $stringency $fitmode

####################### Move to final directory
mkdir -p $outdir
cp -r $out_folder/* $outdir/

####################### CleanUp
rm -r $out_folder
rm -r $tmp_folder
