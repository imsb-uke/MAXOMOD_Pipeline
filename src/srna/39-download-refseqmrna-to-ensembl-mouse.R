#!/usr/bin/env Rscript

cache_biomart <- file.path(getwd(), ".cache_biomaRt")
dir.create(cache_biomart, showWarnings = FALSE, recursive = TRUE)
Sys.setenv("BIOMART_CACHE" = cache_biomart)

library("biomaRt")
ensembl<-  useMart("ensembl", dataset="mmusculus_gene_ensembl", host="nov2020.archive.ensembl.org")

#values<- c("NM_001101", "NM_001256799", "NM_000594")
#getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = values, mart= ensembl)

refseqmrna_ensembl <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "mgi_symbol"), mart= ensembl)
dir.create("database/biomaRt", recursive=TRUE, showWarnings=FALSE)
write.csv(refseqmrna_ensembl, "database/biomaRt/refseq_mrna_to_ensembl_gene_id_mouse.csv", row.names = FALSE)
