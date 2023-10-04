#!/bin/bash

set -e # Fail on error 

#### INFO: to run shellscript, use 'source shellscript.sh'
model=$2     # dataset name (x-mouse/ human)
input_dir=$4 # input directory
gtf=$6      # annotation file
out_dir=$8   # output directory
#https://github.com/mrbaseman/parse_yaml/blob/master/README.md


echo outdir: ${out_dir}
echo input: ${input_dir}
echo gtf: ${gtf}
echo model: ${model}
###################### FEMALE ######################

readarray -t mut_file < ${input_dir}/female_quantfiles_mut.txt
readarray -t ctrl_file < ${input_dir}/female_quantfiles_ctrl.txt

readarray -t mut_file_male < ${input_dir}/male_quantfiles_mut.txt
readarray -t ctrl_file_male < ${input_dir}/male_quantfiles_ctrl.txt

### build file system
mkdir -p ${out_dir}/female/

directory=$(realpath ${out_dir}/female/)

mkdir -p ${directory}/PSI ${directory}/AVGLOGTPM ${directory}/DPSI ${directory}/EVENTS ${directory}/TPM



#### get TPM columns of  quantification files from salmon

python3 src/rnaseq/SUPPA/multipleFieldSelection.py -i ${mut_file[@]} -k 1 -f 4 -o ${directory}/mut_quant_female_${model}.txt
python3 src/rnaseq/SUPPA/multipleFieldSelection.py -i ${ctrl_file[@]} -k 1 -f 4 -o ${directory}/ctrl_quant_female_${model}.txt

#### Generate Events

# is there a way to unpack .gz files temporarily?
#annotation file needs to be a .gtf not compressed (like now)
python3 src/rnaseq/SUPPA/suppa.py generateEvents -i ${gtf} -o ${directory}/Events_${model} -f ioe -e SE SS MX RI FL


#### Calculate PSI per event
# mutant
python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A3_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_A3_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A5_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_A5_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AF_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_AF_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AL_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_AL_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_MX_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_MX_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_RI_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_RI_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_SE_strict.ioe -e ${directory}/mut_quant_female_${model}.txt -o ${directory}/PSI_${model}_SE_mut



# control
python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A3_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_A3_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A5_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_A5_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AF_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_AF_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AL_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_AL_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_MX_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_MX_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_RI_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_RI_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_SE_strict.ioe -e ${directory}/ctrl_quant_female_${model}.txt -o ${directory}/PSI_${model}_SE_ctrl



# get header of event 1, then append all psi values of all events
# put all .ioe events in the same file
head -1 ${directory}/*A3_strict.ioe > ${directory}/all_events.ioe 
tail -q -n +2 ${directory}/*.ioe >> ${directory}/all_events.ioe

head -1 ${directory}/*A3_mut.psi > ${directory}/all_events_mut_together.psi 
tail -q -n +2 ${directory}/*_mut.psi >> ${directory}/all_events_mut_together.psi

head -1 ${directory}/*A3_ctrl.psi > ${directory}/all_events_ctrl_together.psi 
tail -q -n +2 ${directory}/*_ctrl.psi >> ${directory}/all_events_ctrl_together.psi

#### differential splicing analysis for all events
python3 src/rnaseq/SUPPA/suppa.py diffSplice -m empirical -i ${directory}/all_events.ioe -p ${directory}/all_events_mut_together.psi ${directory}/all_events_ctrl_together.psi -e ${directory}/mut_quant_female_${model}.txt ${directory}/ctrl_quant_female_${model}.txt -l 0.05 -gc -o ${directory}/result_allEvents_together --save_tpm_events

mv ${directory}/*.psi ${directory}/PSI
mv ${directory}/*.ioe ${directory}/EVENTS
mv ${directory}/*.gtf ${directory}/EVENTS
mv ${directory}/*avglogtpm.tab ${directory}/AVGLOGTPM
mv ${directory}/*.dpsi* ${directory}/DPSI
mv ${directory}/*.psivec ${directory}/DPSI
mv ${directory}/*.txt ${directory}/TPM

# options: 
#-p = PSI <mut.psi> <ctrl.psi>
#-e = expression/quantification files <mut.txt> <ctrl.txt>
#-l = lower bound (p-value)
#-gc = multiple testing within a gene




###################### MALE

### build file system

mkdir -p ${out_dir}/male/

directory=$(realpath ${out_dir}/male/)


mkdir -p ${directory}/PSI ${directory}/AVGLOGTPM ${directory}/DPSI ${directory}/EVENTS ${directory}/TPM

#### get TPM columns of  quantification files from salmon
python3 src/rnaseq/SUPPA/multipleFieldSelection.py -i ${mut_file_male[@]} -k 1 -f 4 -o ${directory}/mut_quant_male_${model}.txt
python3 src/rnaseq/SUPPA/multipleFieldSelection.py -i ${ctrl_file_male[@]} -k 1 -f 4 -o ${directory}/ctrl_quant_male_${model}.txt

#### Generate Events
python3 src/rnaseq/SUPPA/suppa.py generateEvents -i ${gtf} -o ${directory}/Events_${model} -f ioe -e SE SS MX RI FL


#### Calculate PSI per event
# mutant
python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A3_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_A3_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A5_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_A5_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AF_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_AF_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AL_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_AL_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_MX_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_MX_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_RI_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_RI_mut

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_SE_strict.ioe -e ${directory}/mut_quant_male_${model}.txt -o ${directory}/PSI_${model}_SE_mut



# control
python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A3_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_A3_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_A5_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_A5_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AF_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_AF_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_AL_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_AL_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_MX_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_MX_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_RI_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_RI_ctrl

python3 src/rnaseq/SUPPA/suppa.py psiPerEvent -i ${directory}/Events_${model}_SE_strict.ioe -e ${directory}/ctrl_quant_male_${model}.txt -o ${directory}/PSI_${model}_SE_ctrl



# get header of event 1, then append all psi values of all events
# put all .ioe events in the same file
head -1 ${directory}/*A3_strict.ioe > ${directory}/all_events.ioe 
tail -q -n +2 ${directory}/*.ioe >> ${directory}/all_events.ioe

head -1 ${directory}/*A3_mut.psi > ${directory}/all_events_mut_together.psi 
tail -q -n +2 ${directory}/*_mut.psi >> ${directory}/all_events_mut_together.psi

head -1 ${directory}/*A3_ctrl.psi > ${directory}/all_events_ctrl_together.psi 
tail -q -n +2 ${directory}/*_ctrl.psi >> ${directory}/all_events_ctrl_together.psi

#### differential splicing analysis for all events
python3 src/rnaseq/SUPPA/suppa.py diffSplice -m empirical -i ${directory}/all_events.ioe -p ${directory}/all_events_mut_together.psi ${directory}/all_events_ctrl_together.psi -e ${directory}/mut_quant_male_${model}.txt ${directory}/ctrl_quant_male_${model}.txt -l 0.05 -gc -o ${directory}/result_allEvents_together --save_tpm_events


mv ${directory}/*.psi ${directory}/PSI
mv ${directory}/*.ioe ${directory}/EVENTS
mv ${directory}/*.gtf ${directory}/EVENTS
mv ${directory}/*avglogtpm.tab ${directory}/AVGLOGTPM
mv ${directory}/*.dpsi* ${directory}/DPSI
mv ${directory}/*.psivec ${directory}/DPSI
mv ${directory}/*.txt ${directory}/TPM
