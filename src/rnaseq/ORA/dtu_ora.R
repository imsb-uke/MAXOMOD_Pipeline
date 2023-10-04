#!/usr/bin/env Rscript

library('clusterProfiler')
library('ggplot2')
library(tibble)


############ PREPARATIONS ############
source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/memoise_decorator.R"))
source(here::here("src/shared/R-helpers/enrichGO_memoise.R"))
settings <- get_settings(omic = "splicing_analysis", stage = "ora_dtu", dataset_name = NULL) 

#### Create output directory
if(!dir.exists(settings[["output_directory"]])){
  dir.create(settings[["output_directory"]], recursive = TRUE)
  }
#### params
out_dir <- settings[["output_directory"]]
pval <- settings[['pvalue']]
qval <- settings[['qvalue']]
organism <- settings[['organism']]

  
#### read input files
dex_sign_results_fem <- readr::read_csv(paste0(settings[['inputdex']],"dexseq_FEMALE_padj.csv"), col_names=c("geneID","txID","gene","transcript"), skip=1)
drim_sign_results_fem <- readr::read_csv(paste0(settings[['inputdrim']],"drimseq_FEMALE_padj.csv"), col_names=c("geneID","txID","gene","transcript"), skip=1)
dex_sign_results_male <- readr::read_csv(paste0(settings[['inputdex']],"dexseq_MALE_padj.csv"), col_names=c("geneID","txID","gene","transcript"), skip=1)
drim_sign_results_male <- readr::read_csv(paste0(settings[['inputdrim']],"drimseq_MALE_padj.csv"), col_names=c("geneID","txID","gene","transcript"), skip=1)

######## FEMALE
# significant transcripts out of significant genes
# get significant transcripts out of significant genes
dex_sign_transcripts_fem <- dex_sign_results_fem[dex_sign_results_fem$transcript < 0.05,] 
drim_sign_transcripts_fem <- drim_sign_results_fem[drim_sign_results_fem$transcript < 0.05,]

#union of significant transcripts
dex_IDs_fem <- dplyr::select(dex_sign_transcripts_fem, geneID, txID)
drim_IDs_fem<- dplyr::select(drim_sign_transcripts_fem, geneID, txID)
union_dex_drim_sign_transcripts_fem <- merge(dex_IDs_fem, drim_IDs_fem, all=TRUE)

######### MALE
# significant transcripts out of significant genes
# get significant transcripts out of significant genes
dex_sign_transcripts_male <- dex_sign_results_male[dex_sign_results_male$transcript < 0.05,] 
drim_sign_transcripts_male <- drim_sign_results_male[drim_sign_results_male$transcript < 0.05,]
#union of significant transcripts
dex_IDs_male <- dplyr::select(dex_sign_transcripts_male, geneID, txID)
drim_IDs_male<- dplyr::select(drim_sign_transcripts_male, geneID, txID)
union_dex_drim_sign_transcripts_male <- merge(dex_IDs_male, drim_IDs_male, all=TRUE)



############ ORA ############
gene_union_fem <- unique(union_dex_drim_sign_transcripts_fem$geneID)
gene_union_male<-unique(union_dex_drim_sign_transcripts_male$geneID)
   
if(organism == 'mmusculus'){
	orgdb <- "org.Mm.eg.db"
}else{ orgdb <- "org.Hs.eg.db"}
   
############## female
ORA_fem <- enrichGO_memoise(
            cache_dirname = ".cache_dtu_gene_enrichment/",
            gene = gene_union_fem,
            OrgDb = orgdb,
            keyType = 'ENSEMBL',
            ont = "ALL",
            pAdjustMethod = "BH",
            readable = TRUE,
            pvalueCutoff  = pval,
            qvalueCutoff  = qval
    )


   
############### male
ORA_male <- enrichGO_memoise(
            cache_dirname = ".cache_dtu_gene_enrichment/",
            gene = gene_union_male,
            OrgDb = orgdb,
            keyType = 'ENSEMBL',
            ont = "ALL",
            pAdjustMethod = "BH",
            readable = TRUE,
            pvalueCutoff  = pval,
            qvalueCutoff  = qval
    )

##### save female results
if(nrow(ORA_fem) > 0){
write.csv(ORA_fem,paste0(out_dir,'/ORA_FEM_pval',pval, '_qval',qval,'.csv'), row.names=FALSE)

# calculate plot height depending of number of overrepresented pathways
if(nrow(ORA_fem) > 7 && nrow(ORA_fem) <= 30){fem_heigth <- 3+(nrow(ORA_fem)*0.2)}else if(nrow(ORA_fem) > 30){fem_heigth <- 9.2}else{fem_heigth <- 3}

ORA_female <- setReadable(ORA_fem, orgdb, "ENSEMBL")
if(!any(is.na(ORA_fem$qvalue))){
e1 <- enrichplot::cnetplot(ORA_female, label = 'ALL', cex_label_gene = 0.9, showCategory = 6, cex_label_category = 0.9) + theme( legend.text=element_text(size=14), legend.title = element_text(size=14))
ggsave(filename=paste0(out_dir,"cnetplot_female.pdf"), plot=e1, width=9, height=10, units="in")
    
b1 <- barplot(ORA_fem, showCategory=30, order=TRUE, font.size = 15) + theme(plot.title = element_text(size = 19), legend.text=element_text(size=14), legend.title = element_text(size=14))
print(ggsave(filename=paste0(out_dir,"barplot_female_.pdf"), plot=b1, scale=2, width=4, height=fem_heigth, units="in"))}
}


##### save male results
if(nrow(ORA_male) > 0){
write.csv(ORA_male,paste0(out_dir,'/ORA_MALE_pval',pval, '_qval',qval,'.csv'), row.names=FALSE)
    
# calculate plot height depending of number of overrepresented pathways
    
if(nrow(ORA_male) > 7 && nrow(ORA_male) <= 30){male_heigth <- 3+(nrow(ORA_male)*0.2)}else if(nrow(ORA_male) > 30){male_heigth <- 9.2}else{male_heigth <- 3}
    
    
ORA_male <- setReadable(ORA_male,orgdb, "ENSEMBL")
if(!any(is.na(ORA_male$qvalue))){
e2 <- enrichplot::cnetplot(ORA_male, label = 'ALL', cex_label_gene = 0.9, showCategory = 6, cex_label_category = 0.9) + theme( legend.text=element_text(size=14), legend.title = element_text(size=14))
ggsave(filename=paste0(out_dir,"cnetplot_male.pdf"), plot=e2, width=9, height=10, units="in")

b2 <- barplot(ORA_male, showCategory=30, order=TRUE, font.size = 15)+ theme(plot.title = element_text(size = 19), legend.text=element_text(size=14), legend.title = element_text(size=14))
ggsave(filename=paste0(out_dir,"barplot_male.pdf"), plot=b2,  scale=2, width=4, height=male_heigth, units="in")}
}