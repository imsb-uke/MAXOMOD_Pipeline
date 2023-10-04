#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/gseGO_memoise.R"))
source(here::here("src/shared/R-helpers/gseKEGG_memoise.R"))


suppressPackageStartupMessages({
  library("AnnotationDbi")
  library("dplyr")
  library("purrr")
  library("ggplot2")
  library("enrichplot")
  library("stringr")
})

plot_height <- 30

settings <- get_settings(omic = "phosphoproteomics", stage = "gsea_enrichment", dataset_name = NULL)

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

message("Gene Enrichment analysis...")
if (settings[["organism"]] == "mmusculus") {
  org_db <- "org.Mm.eg.db"
  org_kegg <- "mmu"
} else if (settings[["organism"]] == "hsapiens") {
  org_db <- "org.Hs.eg.db"
  org_kegg <- "hsa"
} else {
  stop("Not implemented error")
}
library(org_db, character.only = TRUE)

de_results <- readRDS(settings[["de_results_rds"]])

message("Reading protein to gene mapping...")
prot_to <- readr::read_tsv(settings[["uniprot_to_ensembl"]], col_names = FALSE)
prot_to_gene <- prot_to[,c("X2", "X19")]
colnames(prot_to_gene) <- c("UniProtName","gene_id")
prot_to_gene_vec <- tibble::deframe(prot_to_gene)
prot_to_gene_vec <- gsub(";.*", "", prot_to_gene_vec)
prot_to_gene_vec <- prot_to_gene_vec[!is.na(prot_to_gene_vec)]

protein_to_gene_id <- function(x, remove_missing = TRUE) {
  out <- prot_to_gene_vec[x]
  names(out) <- x
  if (isTRUE(remove_missing)) {
    out <- out[!is.na(out)]
  }
  out
}

prot_to_entrez <- prot_to[,c("X2", "X3")]
colnames(prot_to_entrez) <- c("UniProtName","gene_id")
prot_to_entrez_vec <- tibble::deframe(prot_to_entrez)
prot_to_entrez_vec <- gsub(";.*", "", prot_to_entrez_vec)
prot_to_entrez_vec <- prot_to_entrez_vec[!is.na(prot_to_entrez_vec)]

protein_to_entrez <- function(x, remove_missing = TRUE) {
  out <- prot_to_entrez_vec[x]
  names(out) <- x
  if (isTRUE(remove_missing)) {
    out <- out[!is.na(out)]
  }
  out
}



message("Starting GSEA calculations...")
gsea_results_dir <- file.path(settings[["output_directory"]], "gsea_results")
dir.create(gsea_results_dir, showWarnings = FALSE, recursive = TRUE)
gsea_results <- purrr::imap(de_results, function(analysis, analysis_name) {
  message(sprintf(" - %s", analysis_name))
  purrr::imap(analysis, function(result, result_name) {
  message(sprintf("  * %s", result_name))
    result <- result %>% 
      mutate(UniProtName=stringr::str_extract(UniProtName,"([A-Z1-9]+_(MOUSE|HUMAN))")) %>%
      filter(!duplicated(UniProtName))
    result$gene <- protein_to_gene_id(result$UniProtName, remove_missing = FALSE)
    result_with_genes <- dplyr::filter(result, !is.na(gene))
    geneList <- result_with_genes$log2FoldChange
    names(geneList) <- result_with_genes$gene
    geneList <- sort(geneList, decreasing = TRUE)
    gsea_results <- c("BP", "MF", "CC") %>%
      purrr::set_names() %>%
      purrr::map(function(ont) {
        message(sprintf("     + %s", ont))
        res <- gseGO_memoise(
          cache_dirname = file.path(settings[["cache"]], "_gseGO_cache"),
          geneList = geneList,
          ont = ont,
          OrgDb = AnnotationDbi::get(org_db),
          keyType = "ENSEMBL",
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500,
          eps = 1e-10,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH"
        )
        if (nrow(as.data.frame(res)) > 0) {
          dp <- enrichplot::dotplot(res, showCategory=30, split=".sign") +
            ggplot2::labs(
              title = sprintf("GO: %s", ont),
              subtitle = sprintf("%s - %s", analysis_name, result_name)
            ) +
            ggplot2::facet_grid( . ~ .sign)
          ggplot2::ggsave(
            dp,
            filename = file.path(
              gsea_results_dir,
              sprintf("%s-%s-%s_dotplot.png", analysis_name, result_name, ont)
            ),
            width = 20,
            height = plot_height,
            units="cm"
          )
        }
        res
      })
    message("     + KEGG")
    result$gene <- protein_to_entrez(result$UniProtName, remove_missing = FALSE)
    result_with_genes <- dplyr::filter(result, !is.na(gene))
    geneList <- result_with_genes$log2FoldChange
    names(geneList) <- result_with_genes$gene
    geneList <- sort(geneList, decreasing = TRUE)
    geneList <- geneList[!duplicated(names(geneList))]
    stopifnot(sum(duplicated(names(geneList))) == 0)
    gsea_results$KEGG <- gseKEGG_memoise(
      cache_dirname = settings[["cache"]],
      geneList = geneList,
      organism = org_kegg,
      keyType = "ncbi-geneid",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )

    if (!is.null(gsea_results$KEGG)&&nrow(gsea_results$KEGG) > 0) {
      dp <- enrichplot::dotplot(gsea_results$KEGG, showCategory=30, split=".sign") +
        ggplot2::labs(
          title = "KEGG",
          subtitle = sprintf("%s - %s", analysis_name, result_name)
        ) +
        ggplot2::facet_grid( . ~ .sign)
      ggplot2::ggsave(
        dp,
        filename = file.path(
          gsea_results_dir,
          sprintf("%s-%s-%s_dotplot.png", analysis_name, result_name, "kegg")
        ),
        width = 20,
        height = plot_height,
        units="cm"
      )
    }
    purrr::imap(gsea_results, function(result, ont) {
      write.csv(
        x = result,
        file = file.path(gsea_results_dir, sprintf("%s-%s-%s_gsea_results.csv", analysis_name, result_name, ont)),
      )
    })
    saveRDS(
      gsea_results,
      file.path(gsea_results_dir, sprintf("%s-%s_gsea_results.rds", analysis_name, result_name))
    )

    gsea_results
  })
})

message("DONE")
