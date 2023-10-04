#!/usr/bin/env Rscript

cache_biomart <- file.path(getwd(), ".cache_biomaRt")
dir.create(cache_biomart, showWarnings = FALSE, recursive = TRUE)
Sys.setenv("BIOMART_CACHE" = cache_biomart)


library("biomaRt")
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl", host="nov2020.archive.ensembl.org")


ensembl_symbol <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart= ensembl)
ensembl_symbol <- dplyr::filter(ensembl_symbol, hgnc_symbol != "")
ensembl_symbol <- ensembl_symbol[!duplicated(ensembl_symbol$ensembl_gene_id),]

dir.create("database/biomaRt", recursive=TRUE, showWarnings=FALSE)
write.csv(ensembl_symbol, "database/biomaRt/ensembl_to_symbol_human.csv", row.names = FALSE)
saveRDS(ensembl_symbol, file = "database/biomaRt/ensembl_to_symbol_human.rds")

