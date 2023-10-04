#!/usr/bin/env Rscript


# multiple hypothesis testing for DTU
StageR_func <- function(pScreen, pConfirmation, tx2gene, pScreenAdjusted){
  stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, tx2gene = tx2gene,
                        pScreenAdjusted = pScreenAdjusted)
  stageRObj <- stageWiseAdjustment(stageRObj, alpha = 0.05, method = "dtu") 
  
  suppressWarnings({
    dtu.p_adj_gene_values <- getAdjustedPValues(stageRObj, order=FALSE,
                                                onlySignificantGenes=TRUE)
  })
  return(dtu.p_adj_gene_values)  
}
