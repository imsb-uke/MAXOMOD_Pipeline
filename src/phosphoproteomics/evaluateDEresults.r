#!/usr/bin/env Rscript

library(argparse)
library(yaml)
library(tidyverse)
library(ggplot2)
library(readxl)
library(cowplot)

#Script Parameters
COLORLEGENDLABEL <- "Measurement ratio to NA"
SIZELEGENDLABEL <- "Number of Measurements"
SCALECOLORGRADIENTLOW <- "Red"
SCALECOLORGRADIENTHIGH <- "#fee0d2" #light red



create_NAratio_volcano<- function(rawInput,DEinput){
    MEratios <- data.frame(MEratios=1-(rowSums(rawInput==0)/length(colnames(rawInput))),
                          MEnumber=length(colnames(rawInput))-rowSums(rawInput==0)) %>%
    mutate(PhosphositeName=rownames(rawInput)) %>%
    inner_join(DEinput, by=c("PhosphositeName"="UniProtName"))

    print(head(MEratios))
    print(colnames(rawInput))

    ggplot(MEratios,aes(x=log2FoldChange,y=-log10(padj),size=MEnumber,color=MEratios)) +
        geom_point() +
        ggtitle(paste("Total number of samples:",length(colnames(rawInput)))) +
        labs(color=COLORLEGENDLABEL,
            size=SIZELEGENDLABEL) +
        scale_color_gradient(low = SCALECOLORGRADIENTLOW, high = SCALECOLORGRADIENTHIGH)
}

get_args <- function() {

    # create parser object
    parser <- ArgumentParser()
    parser$add_argument("--stage-name", type = "character", default = NULL, help = "stage key from params.yaml to use")


    args <- parser$parse_args()
    args
}

# MAIN


## Read params

args = get_args()

stage_name = args$stage_name[1] # selects only first argument for stage_name
print(stage_name)

params = yaml::read_yaml("params.yaml")[["phosphoproteomics"]][["evaluateDEresults"]][["stages"]][[stage_name]]


rawInput<-read_csv(params[["rawInput"]]) %>%
    column_to_rownames(var="PhosphositeName")
print(head(rawInput))

rawAnno <- read_csv(params[["rawAnnotation"]])
print(head(rawAnno))

if(params[["sex"]] %in% unique(rawAnno$Sex)){
    rawAnno <- rawAnno %>%
        filter(Sex==params[["sex"]])
    
    rawInput <- rawInput[,c(rawAnno$SampleID)]
        
} else {
    stop(paste("Not a valid sex! Choose from:", unique(rawAnno$Sex)))
}    

print(rawAnno)

for (group in c("hamburg","strasbourg")){
#LOAD DE input data and get the data from different sources into similar shape
if(str_detect(params[[paste("DEinput",group,sep="_")]],".csv")){
    DEinput<-read_csv(params[[paste("DEinput",group,sep="_")]])
    print(DEinput)
} else if(str_detect(params[[paste("DEinput",group,sep="_")]],".xlsx")){
    if (params[["sex"]]=="female"){
        DEinput<-read_excel(file.path(params[[paste("DEinput",group,sep="_")]]),skip=1) %>%
            rename(log2FoldChange=`logFC (TG F_vs_WT F)`,
              padj=`P_Value (TG F_vs_WT F)`)
        } else if (params[["sex"]]=="male"){
        DEinput<-read_excel(file.path(params[[paste("DEinput",group,sep="_")]]),skip=1) %>%
            rename(log2FoldChange=`logFC (TG M_vs_WT M)`,
              padj=`P_Value (TG M_vs_WT M)`)
    } else {
        stop(paste("Not a valid sex!",params[["sex"]],"Must be of 'female' or 'male'."))
    }
    
    
    DEinput$UniProtName <- paste0(DEinput$Position,
                         DEinput$`Amino acid`,
                         rep('_',nrow(DEinput)),
                         str_extract(DEinput$Protein,"([A-Z0-9]*_[A-Z]*)"),
                         DEinput$Multiplicity)

    print(head(DEinput))
} else {
    stop("Did not recognize file ending! Is it a valid file(.csv, .xlsx)?")
}

#Create NAratio volcano with all samples
create_NAratio_volcano(rawInput,DEinput)
dir.create(params[["outdir"]], recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(params[["outdir"]],paste0(stage_name,'_',group,".pdf")))

#create NAratio plots for all conditions and separated by condition
ggplt_list <- list()
conditions <- unique(rawAnno$Condition)
for (condition in conditions){
    filteredAnno <- rawAnno %>%
        filter(Condition==condition)
    
    print(filteredAnno)

    ggplt_list[[condition]] <- create_NAratio_volcano(rawInput[,c(filteredAnno$SampleID)],
                                                      DEinput)
}

ggsave(file.path(params[["outdir"]],paste0(stage_name,'_',group,"_byCondition.pdf")),
      plot=cowplot::plot_grid(plotlist = ggplt_list,
                             labels=conditions,
                             nrow=2))

}

#CREATE DOTPLOT ON RAW INTENSITIES FOR SOME EXAMPLES
if (!is.null(params[["examples"]])){
examples <- params[["examples"]] 
print(examples)
#print(head(t(rawInput)))
#print(class(t(rawInput)))
#print(pivot_longer(as.data.frame(t(rawInput)),all_of(examples),names_to="Phosphosite",values_to="Intensity"))
ggplot(pivot_longer(as.data.frame(t(rawInput)),all_of(examples),names_to="Phosphosite",values_to="Intensity"),
       aes(x=Intensity))+
    geom_dotplot()+
    scale_x_log10()+
       facet_wrap(~ Phosphosite)
ggsave(file.path(params[["outdir"]],paste0(stage_name,"_dotplotRawIntensitiesOfExamples.pdf")))

    
boxplot_data <- rownames_to_column(as.data.frame(t(rawInput))) %>%
    inner_join(rawAnno, by=c("rowname"="SampleID"))
#print(head(boxplot_data))
#print(str(boxplot_data))
    
BOXPLOT_JITTER_SIZE=0.8
BOXPLOT_ALPHA=0.5
ggplot(pivot_longer(as.data.frame(t(rawInput)),all_of(examples),names_to="Phosphosite",values_to="Intensity"),
       aes(x=1,y=Intensity))+
    geom_boxplot()+
    geom_jitter(size=BOXPLOT_JITTER_SIZE, alpha=BOXPLOT_ALPHA) +
    theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))+
#    scale_y_log10()+
    facet_wrap(~ Phosphosite, scales="free")
ggsave(file.path(params[["outdir"]],paste0(stage_name,"_boxplotRawIntensitiesOfExamplesFaceted.pdf")))

ggplot(pivot_longer(boxplot_data,all_of(examples),names_to="Phosphosite",values_to="Intensity"),
       aes(x=Condition,y=Intensity))+
    geom_boxplot()+
    geom_jitter(size=BOXPLOT_JITTER_SIZE, alpha=BOXPLOT_ALPHA) +
    theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))+
#    scale_y_log10()+
    facet_wrap(~ Phosphosite, scales="free")
ggsave(file.path(params[["outdir"]],paste0(stage_name,"_boxplotRawIntensitiesOfExamplesFacetedbyCondition.pdf")))
    
ggplot(pivot_longer(as.data.frame(t(rawInput)),all_of(examples),names_to="Phosphosite",values_to="Intensity"),
       aes(x=Phosphosite,y=Intensity))+
    geom_boxplot()+
    geom_jitter(size=BOXPLOT_JITTER_SIZE, alpha=BOXPLOT_ALPHA) +
    theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))#+
#    scale_y_log10()+
#    facet_wrap(~ Phosphosite, scales="free")
ggsave(file.path(params[["outdir"]],paste0(stage_name,"_boxplotRawIntensitiesOfExamples.pdf")))
}

#CREATE DOTPLOT OF SAMPLES AFTER IMPUTATION PER HIT

if(str_detect(params[["downstreamInput"]],".csv")){
    downstreamInput<-t(read.csv(params[["downstreamInput"]], row.names=1))
    
} else if(str_detect(params[["downstreamInput"]],".xlsx")){
    downstreamInput<-read_excel(file.path(params[["downstreamInput"]]))
} else {
    stop("Did not recognize file ending! Is it a valid file (.csv, .xlsx)?")
}

print(head(downstreamInput))

#CREATE table with cols :phosphoprotein,
#number of ctrl samples imputed,
#number of mutated samples imputed,
#number of ctrl samples available,
#number of mutated samples available,
#log2FC strasbourg,
#log2FC hamburg,
#padj for strasbourg,
#padj for hamburg

DE_hamburg <- read_csv(params[[paste("DEinput_hamburg")]])
if (params[["sex"]]=="female"){
    DE_strasbourg <- read_excel(file.path(params[[paste("DEinput",group,sep="_")]]),skip=1) %>%
        rename(log2FoldChange=`logFC (TG F_vs_WT F)`,
              padj=`P_Value (TG F_vs_WT F)`,
              isDifferential=`isDifferential (TG F_vs_WT F)`)
        } else if (params[["sex"]]=="male"){
    DE_strasbourg <- read_excel(file.path(params[[paste("DEinput",group,sep="_")]]),skip=1) %>%
        rename(log2FoldChange=`logFC (TG M_vs_WT M)`,
              padj=`P_Value (TG M_vs_WT M)`,
              isDifferential=`isDifferential (TG M_vs_WT M)`)
    } else {
        stop(paste("Not a valid sex!",params[["sex"]],"Must be of 'female' or 'male'."))
    }

DE_strasbourg$UniProtName <- paste0(DE_strasbourg$Position,
                         DE_strasbourg$`Amino acid`,
                         rep('_',nrow(DE_strasbourg)),
                         str_extract(DE_strasbourg$Protein,"([A-Z0-9]*_[A-Z]*)"),
                         DE_strasbourg$Multiplicity)

joined_DE <- DE_hamburg %>%
    full_join(DE_strasbourg, by="UniProtName" , suffix=c("_hamburg","_strasbourg")) %>%
    select(UniProtName, 
           padj_hamburg, 
           padj_strasbourg, 
           log2FoldChange_hamburg, 
           log2FoldChange_strasbourg, 
           isDifferential) %>%
    filter(isDifferential==1)

print(str(joined_DE))
print(head(joined_DE))

joined_DE <- data.frame(numCTRLimputed=rowSums(rawInput[,rawAnno$SampleID[rawAnno$Condition=="ctrl"]]==0),
                        numMUTimputed=rowSums(rawInput[,rawAnno$SampleID[rawAnno$Condition=="mut"]]==0),
                        numCTRL=rep(length(rawAnno$SampleID[rawAnno$Condition=="ctrl"]),times=nrow(rawInput)),
                        numMUT=rep(length(rawAnno$SampleID[rawAnno$Condition=="mut"]),times=nrow(rawInput))) %>%
    #MEratios=1-(rowSums(rawInput==0)/length(colnames(rawInput))),
    #                  MEnumber=length(colnames(rawInput))-rowSums(rawInput==0)) %>%
    mutate(PhosphositeName=rownames(rawInput)) %>%
    inner_join(joined_DE, by=c("PhosphositeName"="UniProtName")) %>%
    select(PhosphositeName, numCTRLimputed, numMUTimputed, numCTRL, numMUT, padj_hamburg, padj_strasbourg,log2FoldChange_hamburg,log2FoldChange_strasbourg, -isDifferential)

print(str(joined_DE))
print(head(joined_DE))

write.csv(joined_DE, file.path(params[["outdir"]],paste0(stage_name,"_DEcomparison.csv")), row.names = FALSE)