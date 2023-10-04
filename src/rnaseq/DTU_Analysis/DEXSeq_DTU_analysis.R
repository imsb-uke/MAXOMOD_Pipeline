#!/usr/bin/env Rscript
#################### Differential Transcript Usage Analysis with DEXSeq
library(DEXSeq)

source(here::here("src/rnaseq/DTU_Analysis/stageR_DTU_significanceTesting.R"))

DEXSeq_method <- function(sample_info, filtered_dm_object, design_full, sex, organism){

if(organism == 'mmusculus'){
strp <- function(x) substr(x,1,18) 
}else{
strp <- function(x) substr(x,1,15)}
  sample.data <- sample_info    
  count.data <- round(as.matrix(counts(filtered_dm_object)[,-c(1:2)])) # round counts and put counts without ID in count.data
  dxd <- DEXSeqDataSet(countData=count.data,
                       sampleData=sample.data,
                       design=~sample + exon + Condition:exon, # differences in exon usage due to changes in the 'condition' 
                       featureID=counts(filtered_dm_object)$feature_id,
                       groupID=counts(filtered_dm_object)$gene_id)
  
  cat(" - Starting DEXSeq analysis - ", file=stdout())
  system.time({
    dxd <- estimateSizeFactors(dxd) #Normalisation?
    dxd <- estimateDispersions(dxd, quiet=TRUE) #estimate variability of the data (Streuung), parameters and kurve
    dxd <- testForDEU(dxd, reducedModel=~sample + exon) #likelihood ratio test for differential transcript usage 
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "Condition")  
  })
    cat(" - Finished DEXSeq analysis - ")

  
  dex_results <- DEXSeqResults(dxd, independentFiltering = FALSE )
  full_results <- dex_results
  
  #### StageR
  #Compute per-gene adjusted p-value (aggregates evidence from multiple tests within 
  # a gene to a single p-value for the gene and then corrects for multiple testing across genes)
  q_values <- perGeneQValue(dex_results) 
  columns <- c("featureID","groupID","pvalue")
  dex_results_short <- as.data.frame(dex_results[,columns])
  
  # pScreen = Vector with p-values and gene IDs
  pScreen <- q_values   
  names(pScreen) <- strp(names(pScreen))   #strip gene IDs to 18 characters
  
  # pConfirmation = p-values with transcript_IDs
  pConfirmation <- matrix(dex_results_short$pvalue, ncol=1)
  dimnames(pConfirmation) <- list(strp(dex_results_short$featureID))
  
  # tx2gene = 2-column data.frame with transcript and gene IDs (stripped)
  tx2gene <- as.data.frame(dex_results_short[,c("featureID", "groupID")])
  tx2gene$groupID <- strp(tx2gene$groupID) #strip group ID to 18 characters
  tx2gene$featureID <- strp(tx2gene$featureID)
  
  
  dex.p_adj_values <- StageR_func(pScreen, pConfirmation, tx2gene, pScreenAdjusted = TRUE)
 
  
  return(dex.p_adj_values)
}
