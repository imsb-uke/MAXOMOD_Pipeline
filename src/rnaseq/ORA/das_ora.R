#!/usr/bin/env Rscript


############ PREPARATIONS ############
library('clusterProfiler')
library('biomaRt') 
library('stringr')
library('ggplot2')
source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/memoise_decorator.R"))
source(here::here("src/shared/R-helpers/enrichGO_memoise.R"))

settings <- get_settings(omic = "splicing_analysis_das", stage = "das_ora", dataset_name = NULL) 

Sys.setenv(BIOMART_CACHE = ".cache_biomart_gene_enrichment/")

# Create output directory
if(!dir.exists(settings[["output_directory"]])){
  dir.create(settings[["output_directory"]], recursive = TRUE)
  }

# params
out_dir <- settings[["output_directory"]]
pval <- settings[['pvalue']]
qval <- settings[['qvalue']]
organism <- settings[['organism']]

cnetplot <- function(res) {
    out <- tryCatch(
        {
           out <- enrichplot::cnetplot(res, label = 'ALL', cex_label_gene = 0.9, showCategory = 6, cex_label_category = 0.9) + theme( legend.text=element_text(size=14), legend.title = element_text(size=14))
        },
        error=function(cond) {
            message(cond)
            return(NULL)
        }
    )    
    return(out)
}

barplot_catch <- function(res) {
    out <- tryCatch(
        {
           out <- barplot(res, showCategory=30, order=TRUE, font.size = 15) + theme(plot.title = element_text(size = 19), legend.text=element_text(size=14), legend.title = element_text(size=14))
        },
        error=function(cond) {
            message(cond)
            return(NULL)
        }
    )    
    return(out)
}

setReadable_catch <- function(res, orgdb, key) {
    out <- tryCatch(
        {
           out <- setReadable(res, orgdb, key)
        },
        error=function(cond) {
            message(cond)
            return(res)
        }
    )    
    return(out)
}


################## FUNCTIONs

# all_genenames for converting Gene-IDs into GeneSymbols
all_genenames <- function(all_significant_events){
    all_genenames <- gsub( '\\..*', '', rownames(all_significant_events))

    all_unique_genenames <- unique(all_genenames)

    if(organism == 'hsapiens'){set <- 'hsapiens_gene_ensembl'} else {set <- 'mmusculus_gene_ensembl'}                          
    mart_ensembl<- biomaRt::useMart("ensembl", dataset=set, host = "https://nov2020.archive.ensembl.org") 
    # needs to be removed/edited sometimes get the biotype and Gene name of significant genes
    all_geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), filters = "ensembl_gene_id", values = all_unique_genenames , mart = mart_ensembl)
    colnames(all_geneNames) <- c('GeneID', 'GeneSymbol','Description')

    shortGeneID_dpsi_all <- all_significant_events
    shortGeneID_dpsi_all$GeneID <- all_genenames
    full_sign_table_all <- merge(shortGeneID_dpsi_all, all_geneNames)
    
    full_sign_table_all$Event <- rownames(shortGeneID_dpsi_all)
    full_sign_table_all$Event <- str_extract(full_sign_table_all$Event, ";[:alpha:]+[0-9]?")
    full_sign_table_all$Event <- sub(";", "",full_sign_table_all$Event)
    ## reorder columns
    full_sign_table_all <- full_sign_table_all[,c(1,6,2,3,4,5)]
    return(full_sign_table_all)
}
  

################## MAIN

names <- c('dPSI', 'pvalue')
######## get FEMALE significant Events and their GeneSymbol
fem_dpsi <- read.table(paste0(settings[['input']],"female/DPSI/result_allEvents_together.dpsi"))
print(head(fem_dpsi))
colnames(fem_dpsi) <- names
fem_significant <- fem_dpsi[(fem_dpsi$pvalue < 0.05),]
fem_genenames <- all_genenames(fem_significant)
write.csv(fem_genenames, file.path(out_dir, "sign_Events_FEMALE.csv"), row.names = FALSE)

######### get MALE significant Events and their GeneSymbol
male_dpsi <- read.table(paste0(settings[['input']],"male/DPSI/result_allEvents_together.dpsi"))
colnames(male_dpsi) <- names
male_significant <- male_dpsi[(male_dpsi$pvalue < 0.05),]
male_genenames <- all_genenames(male_significant)
write.csv(male_genenames, file.path(out_dir, "sign_Events_MALE.csv"), row.names = FALSE)
################ ORA 
   
if(organism == 'mmusculus'){
	orgdb <- "org.Mm.eg.db"
}else{ orgdb <- "org.Hs.eg.db"}
   
## female
ORA_fem <- enrichGO_memoise(
            cache_dirname = ".cache_das_gene_enrichment/",
            gene = fem_genenames$GeneID,
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
            cache_dirname = ".cache_das_gene_enrichment/",
            gene = male_genenames$GeneID,
            OrgDb = orgdb,
            keyType = 'ENSEMBL',
            ont = "ALL",
            pAdjustMethod = "BH",
            readable = TRUE,
            pvalueCutoff  = pval,
            qvalueCutoff  = qval
    )



##### save female results
print('starting female')
if(nrow(ORA_fem) > 0){
write.csv(ORA_fem,paste0(out_dir,'/ORA_FEM_pval',pval, '_qval',qval,'.csv'), row.names=FALSE)

# calculate plot height depending of number of overrepresented pathways
if(nrow(ORA_fem) > 7 && nrow(ORA_fem) <= 30){fem_heigth <- 3+(nrow(ORA_fem)*0.2)}else if(nrow(ORA_fem) > 30){fem_heigth <- 9.2}else{fem_heigth <- 3}

ORA_female <- setReadable_catch(ORA_fem, orgdb, "ENSEMBL")
print('printing female cnetplot')
if(!any(is.na(ORA_fem$qvalue))){
e1 <- cnetplot(ORA_female)

if (!is.null(e1)){    
print(ggsave(filename=paste0(out_dir,"cnetplot_female.pdf"), plot=e1, width=9, height=10, units="in"))
    }
   
print('printing female barplot')

b1 <- barplot_catch(ORA_fem)
if (!is.null(b1)){    
    print(ggsave(filename=paste0(out_dir,"barplot_female_.pdf"), plot=b1, scale=2, width=4, height=fem_heigth, units="in"))
    }
}
}
print('starting male')

##### save male results
if(nrow(ORA_male) > 0){
write.csv(ORA_male,paste0(out_dir,'/ORA_MALE_pval',pval, '_qval',qval,'.csv'), row.names=FALSE)

print("save done")
    
# calculate plot height depending of number of overrepresented pathways
    
if(nrow(ORA_male) > 7 && nrow(ORA_male) <= 30){male_heigth <- 3+(nrow(ORA_male)*0.2)}else if(nrow(ORA_male) > 30){male_heigth <- 9.2}else{male_heigth <- 3}
    

    
    
ORA_male <- setReadable_catch(ORA_male,orgdb, "ENSEMBL")
if(!any(is.na(ORA_male$qvalue))){
e2 <- cnetplot(ORA_male)
if (!is.null(e2)){  
tryCatch(
        {
           print(ggsave(filename=paste0(out_dir,"cnetplot_male.pdf"), plot=e2, width=9, height=10, units="in"))
        },
        error=function(cond) {
            message(cond)
        }
    )    
    }
    
b2 <- barplot_catch(ORA_male)
if (!is.null(b2)){    
    print(ggsave(filename=paste0(out_dir,"barplot_male.pdf"), plot=b2,  scale=2, width=4, height=male_heigth, units="in"))
    }
}
}