#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

suppressPackageStartupMessages({
  library("clusterProfiler")
  library("ggplot2")
  library("org.Hs.eg.db")
  library("stringr")
})

settings <- get_settings(omic = "integration", stage = "cluster_profiler_compare", dataset_name = NULL) # dataset_name="human"

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

data <- read.csv(settings[["inputs"]])
data <- data[!is.na(data$ENTREZ), ]
data$abs_log2FoldChange <- abs(data$log2FoldChange)

data$ENTREZ <- as.factor(data$ENTREZ)
data$REGULATION <- as.factor(data$REGULATION)

data <- data[data$padj <= settings[["pvalue_cutoff"]], ]

data <- data[data$abs_log2FoldChange >= settings[["absfoldchange_cutoff"]], ]

data <- data[data$log2FoldChange >= settings[["pos_foldchange_cutoff"]], ]

data <- data[data$log2FoldChange <= settings[["neg_foldchange_cutoff"]], ]


formula <- as.formula(paste("ENTREZ", settings[["formula"]], sep = " ~ "))

formula_res <- compareCluster(formula, data = data, fun = settings[["function"]], OrgDb = org.Hs.eg.db)

formula_res <- setReadable(formula_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(formula_res, file = file.path(settings[["output_directory"]], "enrichment.csv"))




if (!is.null(settings[["splitby"]])){
    form <- as.formula(paste0("~",settings[["splitby"]][[2]]))
    p <- dotplot(formula_res, label_format = 10000, x=settings[["splitby"]][[1]]) 
    p <- p + facet_grid(form)
    
    } else {
    p <- dotplot(formula_res, label_format = 10000)
    }
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(p, file = file.path(settings[["output_directory"]], "dotplot.png"), width = 15, height = 8)
ggsave(p, file = file.path(settings[["output_directory"]], "dotplot.pdf"), width = 15, height = 8)

p <- cnetplot(formula_res)
ggsave(p, file = file.path(settings[["output_directory"]], "cnetplot.png"), width = 15, height = 15)
ggsave(p, file = file.path(settings[["output_directory"]], "cnetplot.pdf"), width = 15, height = 15)