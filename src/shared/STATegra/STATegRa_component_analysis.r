#This Script is supposed to analyze multiple omics datasets for shared components in order to decide whcih combinations to use later on

#BiocManager::install("STATegRa")

library(STATegRa)
library(yaml)
library(tidyverse)
library(ggplot2)
library(gridExtra)
#library(grid)

print("yeah! Lets start STATegra!")

params <- yaml::read_yaml("params.yaml")[["stategra"]]

### Transcriptomics data
transcriptomics_exp_mat <- read.csv(params[["stages"]][["c9orf72"]][["omics"]][["transcriptomics"]][["expression_file"]],
                           row.names=1)
#print(dim(transcriptomics_exp_mat))
#print(head(transcriptomics_exp_mat))
transcriptomics_anno <- read.csv(params[["stages"]][["c9orf72"]][["omics"]][["transcriptomics"]][["annotation_file"]],
                           row.names=1)
#print(transcriptomics_anno)

# Block1 - gene expression data
transcriptomics <- createOmicsExpressionSet(Data=as.matrix(transcriptomics_exp_mat), pData=transcriptomics_anno)


### Proteomic data
proteomics_exp_mat <- t(read.csv(params[["stages"]][["c9orf72"]][["omics"]][["proteomics"]][["expression_file"]],
                           row.names=1))
#print(dim(proteomics_exp_mat))
#print(head(proteomics_exp_mat))
proteomics_anno <- read.csv(params[["stages"]][["c9orf72"]][["omics"]][["proteomics"]][["annotation_file"]],
                           row.names=1)
#print(proteomics_anno)

# Block2 - protein expression data
proteomics <- createOmicsExpressionSet(Data=as.matrix(proteomics_exp_mat), pData=proteomics_anno)

## Model selection
for (Rmax in 2:6){
    ms <- modelSelection(Input=list(transcriptomics, proteomics), Rmax=Rmax, fac.sel="single%",
                     varthreshold=0.03, center=TRUE, scale=TRUE, 
                     weight=TRUE, plot_common=FALSE, plot_dist=FALSE)
    ggsave(file.path(params[["out_dir"]],paste0("model_Selection_Rmax",Rmax,".png")),
          grid.arrange(ms$common$pssq, ms$common$pratios, ncol=2) #switching plot_common=TRUE gives automaticaly these plot
       )
    print(names(ms$common))
    print(names(ms))
    
    print(ms$dist$numComps)
    print(ms$common$commonComps)
}