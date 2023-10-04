#!/usr/bin/env Rscript

#################### FUNCTIONS
# function to strip after 18 or 15 characters (remove Ensembl Version ID)    

library(tximport)             
library(GenomicFeatures)      # TxDb Database
# library(DBI)   only if transcript-gene mapping should be saved as sqlite file 
library(DRIMSeq)
library(stageR)
library(biomaRt)


source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/rnaseq/DTU_Analysis/DEXSeq_DTU_analysis.R"))
source(here::here("src/rnaseq/DTU_Analysis/DRIMSeq_DTU_analysis.R"))

#################### Preparations
#import setting file 
settings <- get_settings(omic = "splicing_analysis", stage = "dtu", dataset_name = NULL) # for example: dataset_name="human_DRIMSeq" (?)

# Create output directory
if(!dir.exists(settings[["output_directory"]])){
  dir.create(settings[["output_directory"]], recursive = TRUE)
  }
out_dir <- settings[["output_directory"]]
# Method (DEXSeq or DRIMSeq)
method <- settings[["method"]]
num_samples <-   settings[["num_samples"]]
organism <- settings[["organism"]]

#################### import files from settings (yaml file)
# cohort
sample_info <- readr::read_csv(settings[["cohort"]], col_types = readr::cols())
sample_info <- as.data.frame(sample_info)
print(head(sample_info))

sample_info$SampleID <- gsub('-', '.', sample_info$SampleID)
print(head(sample_info))
nf_ids_to_samples <- readr::read_csv(settings[["nfs_to_sample_id"]], col_types = readr::cols())

#quantification files
path_to_quant <- settings[["input_directory"]]
dirs<- grep("*",list.dirs(path_to_quant,recursive=FALSE),value=TRUE)                            
quant_files<- list.files(dirs, "quant.sf", recursive=TRUE, full.names=TRUE)

if(length(quant_files)<num_samples){
  print("error: not enough quant.sf files found.", file=stderr())
  return(1)
}

# assign sample IDs in the same order as the files
names(quant_files) <- nf_ids_to_samples[order(nf_ids_to_samples$nf_id), ]$SampleID

# reorder quant_files 
quant_files <- quant_files[order(names(quant_files))]

# factor Data (categorize the data and store it as levels)
sample_info$Condition <- factor(sample_info$Condition)
sample_info$Condition<-relevel(sample_info$Condition, "ctrl") 
sample_info$Sex <- factor(sample_info$Sex)

sample_info_female <- sample_info[grep("female", sample_info$Sex),]
sample_info_male <- sample_info[which(sample_info$Sex != "female"),]


#################### tximport (import quantification files into R and scale)
# idea: countsFromAbundance = 'dtuscaledTPM'
tximport_counts <- tximport::tximport(quant_files, type="salmon", txOut=TRUE,countsFromAbundance ="scaledTPM")

#remove all transcripts that are never expressed
cts <- tximport_counts$counts
counts <- cts[rowSums(cts) > 0,]


#################### transcript to gene mapping (TxDb)
gene_annotation_file <-settings[["gtf"]]
txdb_database <- GenomicFeatures::makeTxDbFromGFF(gene_annotation_file, format = "gtf")
print(head(sample_info_female))
print(head(sample_info_male))
txdb_dataframe <- AnnotationDbi::select(txdb_database, keys(txdb_database, "GENEID"), "TXNAME", "GENEID")

## control, if all transcripts are in TxDb Database
if(!(all(rownames(counts) %in% txdb_dataframe$TXNAME))){
  print("error: Not all transcripts in TxDb Database. Please check ensembl ID format.", file=stderr())
  return(1)
}

# select all txdb transcripts, that are in sample data
txdb_similar <- txdb_dataframe[match(rownames(counts),txdb_dataframe$TXNAME),]

if(!(all(rownames(counts) == txdb_similar$TXNAME))){
  print("error: Transcripts of the counts and the TxDb database are different.", file=stderr())
  return(1)
}


# build dataframe of counts and txdb infos (Gene IDs)
counts_and_txdb <- data.frame(gene_id = txdb_similar$GENEID, feature_id = txdb_similar$TXNAME, counts)


#################### Filter transcripts with DRIMSeq filter (dmfilter)
print(head(counts_and_txdb))

colnames(sample_info_female)[colnames(sample_info_female) == "SampleID"] <- "sample_id" 

### FEMALE
      
dm_data_object_female <- DRIMSeq::dmDSdata(counts=counts_and_txdb, samples = sample_info_female) 


#number of mutants
if(organism == 'mmusculus'){
num_mut_fem <- nrow(sample_info_female[grep("mut", sample_info_female$Condition),])
num_mut_male <- nrow(sample_info_male[grep("mut", sample_info_male$Condition),])
}else{
num_mut_fem <- nrow(sample_info_female[grep("als", sample_info_female$Condition),])
num_mut_male <- nrow(sample_info_male[grep("als", sample_info_male$Condition),])}

# number of controls
num_ctrl_fem <- nrow(sample_info_female[grep("ctrl", sample_info_female$Condition),])
num_ctrl_male <- nrow(sample_info_male[grep("ctrl", sample_info_male$Condition),])

n_fem <- (num_samples - (num_mut_male+num_ctrl_male))
print(n_fem) # total number of female samples
n_male <- (num_samples - (num_mut_fem+num_ctrl_fem))
print(cat('\n n male:', n_male))
n.small_fem <- min(c(num_mut_fem,num_ctrl_fem))    # sample size of smallest group (ctrl/mut)
n.small_male <- min(c(num_mut_male,num_ctrl_male))    # sample size of smallest group (ctrl/mut)

filtered_dm_object_female <- DRIMSeq::dmFilter(dm_data_object_female,
min_samps_feature_expr=n.small_fem,min_feature_expr=10,       min_samps_feature_prop=n.small_fem, min_feature_prop=0.1,   min_samps_gene_expr=n_fem, min_gene_expr=10)


### MALE 
colnames(sample_info_male)[colnames(sample_info_male) == "SampleID"] <- "sample_id" 
dm_data_object_male <- DRIMSeq::dmDSdata(counts=counts_and_txdb, samples = sample_info_male) 

filtered_dm_object_male <- DRIMSeq::dmFilter(dm_data_object_male,
                                    min_samps_feature_expr=n.small_male, min_feature_expr=10,
                                    min_samps_feature_prop=n.small_male, min_feature_prop=0.1,
                                    min_samps_gene_expr=n_male, min_gene_expr=10)

print(sample_info_female)
print(sample_info_male)
# FEMALE
all_tested_female <- data.frame(gene_ID <- counts(filtered_dm_object_female)$gene_id, 
                         tx_ID <- counts(filtered_dm_object_female)$feature_id)
colnames(all_tested_female) <- c("gene_ID", "tx_ID")
write.csv(all_tested_female, paste0(out_dir, "/tested_genes_and_features-FEMALE.csv"), row.names= FALSE)

# MALE
all_tested_male <- data.frame(gene_ID <- counts(filtered_dm_object_male)$gene_id, 
                         tx_ID <- counts(filtered_dm_object_male)$feature_id)
colnames(all_tested_male) <- c("gene_ID", "tx_ID")
write.csv(all_tested_male, paste0(out_dir, "tested_genes_and_features-MALE.csv"), row.names= FALSE)


#### Design-matrix
design_full_female <- model.matrix(~Condition, data=DRIMSeq::samples(filtered_dm_object_female))
   
design_full_male <- model.matrix(~Condition, data=DRIMSeq::samples(filtered_dm_object_male))



print(design_full_female)
print(design_full_male)

#################### DRIM or DEX analysis and StageR multiple hypothesis correction
# result: Dataframe with significant Genes and p-adj values for Genes and Transcripts

   
result_female <- get(method)(sample_info_female, filtered_dm_object_female, design_full_female, "FEMALE", organism)
result_male <- get(method)(sample_info_male, filtered_dm_object_male, design_full_male, "MALE", organism)

sign_result_fem <- result_female[result_female$transcript < 0.05,] 

sign_result_male <- result_male[result_male$transcript < 0.05,] 
   
if(method == 'DEXSeq_method'){
  write.csv(result_female, paste0(out_dir,'dexseq_FEMALE_padj.csv'),row.names = FALSE)
  write.csv(result_male, paste0(out_dir,'dexseq_MALE_padj.csv'),row.names = FALSE)
 write.csv(sign_result_male, paste0(out_dir,'sign_dexseq_MALE_padj0.05.csv'),row.names = FALSE)
     write.csv(sign_result_fem, paste0(out_dir,'sign_dexseq_FEMALE_padj0.05.csv'),row.names = FALSE)
    
	  print("DTU analysis finished", file = stdout())

}

if(method == 'DRIMSeq_method'){
  write.csv(result_female, paste0(out_dir,'drimseq_FEMALE_padj.csv'),row.names = FALSE)
  write.csv(result_male, paste0(out_dir,'drimseq_MALE_padj.csv'),row.names = FALSE)
    
     write.csv(sign_result_male, paste0(out_dir,'sign_drimseq_MALE_padj0.05.csv'),row.names = FALSE)
     write.csv(sign_result_fem, paste0(out_dir,'sign_drimseq_FEMALE_padj0.05.csv'),row.names = FALSE)
    
	  print("DTU analysis finished", file = stdout())
}





