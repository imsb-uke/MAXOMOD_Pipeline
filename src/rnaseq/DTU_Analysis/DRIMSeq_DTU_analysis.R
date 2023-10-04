#!/usr/bin/env Rscript
#################### Differential Transcript Usage Analysis with DEXSeq
library(DRIMSeq)

source(here::here("src/rnaseq/DTU_Analysis/stageR_DTU_significanceTesting.R"))


DRIMSeq_method <- function(sample_info, filtered_dm_object, design_full, sex, organism){

if(organism == 'mmusculus'){
strp <- function(x) substr(x,1,18) 
}else{
strp <- function(x) substr(x,1,15)}
  #### DRIMSEQ
  set.seed(1) # random number generator, to make analysis reproducible 
  if(organism=='hsapiens'){coef <- 'Conditionals'
    }else{coef<- 'Conditionmut'}
  cat("- Starting DRIMSeq analysis - ")
  system.time({
    filtered_dm_object <- dmPrecision(filtered_dm_object, design=design_full) # Estimate the precision
    filtered_dm_object <- dmFit(filtered_dm_object, design=design_full) # fit regression coefficient
    filtered_dm_object <- dmTest(filtered_dm_object, coef=coef) # null hypothesis testing on the coefficient of interest
  })
  cat(" - Finished DRIMSeq analysis - ")

  results_gene <- DRIMSeq::results(filtered_dm_object)
  results_txt <- DRIMSeq::results(filtered_dm_object, level = "feature")
  
  #remove NA-values if existing and replace them with '1'
  no.na <- function(x) ifelse(is.na(x), 1, x) 
  results_gene$pvalue <- no.na(results_gene$pvalue)
  results_txt$pvalue <- no.na(results_txt$pvalue)
  
  
  #### StageR
  # p-Screen = Construct a vector of per-gene p-values 
  pScreen <- results_gene$pvalue    # vector with p-values
  names(pScreen) <- strp(results_gene$gene_id)    # assign stripped gene-ids to pvalues
  
  # p-confirmation = p-values with transcript_IDs
  pConfirmation <- matrix(results_txt$pvalue, ncol=1) # p-values with transcript_ids
  rownames(pConfirmation) <- strp(results_txt$feature_id)
  tx2gene <- results_txt[,c("feature_id", "gene_id")]
  
  # strip gene_id to 18 or 15 characters
  tx2gene$gene_id <- strp(tx2gene$gene_id)
  tx2gene$feature_id <- strp(tx2gene$feature_id)
  
  #stageR result (significant Genes)
  drim.padj <- StageR_func(pScreen, pConfirmation, tx2gene, pScreenAdjusted = FALSE)
  
  return(drim.padj)
}


