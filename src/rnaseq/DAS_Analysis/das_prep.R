#!/usr/bin/env Rscript

library(stringr)
library(dplyr)

source(here::here("src/shared/R-helpers/get_settings.R"))

settings <- get_settings(omic = "splicing_analysis_das", stage = "das_prep", dataset_name = NULL) 

if(!dir.exists(settings[['output_directory']])){
    dir.create(settings[['output_directory']], recursive = TRUE)
}
model <- settings[["model"]]
nextflow <- read.csv(settings[["nextflow"]])
annotation <-read.csv(settings[["annotation"]])
quant_file_path <- settings[["input_directory"]]
out_dir <-  settings[['output_directory']]
full_annotation_info <- merge(nextflow, annotation, by='SampleID')
                      
#full_annotation_info$nf_id[which(full_annotation_info$nf_id == 'mut_R8')] <- 'ctrl_R11' #SOD1
                      
nrow(full_annotation_info)
tail(full_annotation_info,50)

female_list <- full_annotation_info[which(full_annotation_info$Sex =='female'),]
male_list <- full_annotation_info[which(full_annotation_info$Sex !='female'),]

   
if (model == 'human-datasets'){
    female_list_mut <- female_list[which(female_list$Condition == 'als'),]
    male_list_mut <- male_list[which(male_list$Condition == 'als'),]

}else{
    female_list_mut <- female_list[which(female_list$Condition == 'mut'),]
    male_list_mut <- male_list[which(male_list$Condition == 'mut'),]
} 

female_list_ctrl <- female_list[which(female_list$Condition == 'ctrl'),]
male_list_ctrl <- male_list[which(male_list$Condition == 'ctrl'),]

### saves all paths to variable
female_quantfiles_mut <- paste0(quant_file_path,female_list_mut$nf_id,'/quant.sf')
female_quantfiles_ctrl <- paste0(quant_file_path,
                                 female_list_ctrl$nf_id, '/quant.sf')

male_quantfiles_mut <- paste0(quant_file_path,male_list_mut$nf_id,'/quant.sf')
male_quantfiles_ctrl <- paste0(quant_file_path,male_list_ctrl$nf_id,'/quant.sf')


### save files 
write.csv(full_annotation_info, paste0(out_dir,'full_annotation.csv'),
          row.names=FALSE)
write.table(female_quantfiles_mut, paste0(out_dir,'female_quantfiles_mut.txt'),
            sep = ' ', row.names=FALSE, col.names =FALSE, quote = FALSE)
write.table(male_quantfiles_mut, paste0(out_dir,'male_quantfiles_mut.txt'),
            sep = ' ', row.names=FALSE, col.names =FALSE, quote = FALSE)
write.table(female_quantfiles_ctrl,paste0(out_dir,'female_quantfiles_ctrl.txt'),
            sep = ' ', row.names=FALSE, col.names =FALSE, quote = FALSE)
write.table(male_quantfiles_ctrl, paste0(out_dir,'male_quantfiles_ctrl.txt'),
            sep = ' ', row.names=FALSE, col.names =FALSE, quote = FALSE)
          
          