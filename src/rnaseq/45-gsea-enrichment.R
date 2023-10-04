#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/remove_ensembl_version.R"))
source(here::here("src/shared/R-helpers/gseGO_memoise.R"))

suppressPackageStartupMessages({
  library("AnnotationDbi")
  library("biomaRt")
  library("clusterProfiler")
  library("dplyr")
  library("ggplot2")
  library("purrr")
  library("readr")
  library("tibble")
})

settings <- get_settings(omic = "rnaseq", stage = "gene_enrichment", dataset_name = NULL) # dataset_name="human")

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

organism_to_annotdbistr <- function(org) {
  if (org == "mmusculus") {
    org_db <- "org.Mm.eg.db"
  } else if (org == "hsapiens") {
    org_db <- "org.Hs.eg.db"
  } else {
    stop("Not implemented error")
  }
  library(org_db, character.only = TRUE)
  org_db
}

org_db <- organism_to_annotdbistr(settings[["organism"]])
orgdb_obj <- AnnotationDbi::get(org_db)

message("Reading DEG outputs")

de_results <- list()
for (deseq2_setting in settings[["gene_enrichment"]]) {
  deseq2_result_dir <- file.path(settings[["deg_dir"]], deseq2_setting[["name"]])
  result_tables <- list()
  for (result in deseq2_setting[["results"]]) {
    result_name <- result[["name"]]
    csv_file <- file.path(deseq2_result_dir, sprintf("%s.csv", result_name))
    de_table <- readr::read_csv(
      file = csv_file,
      col_types = cols(
        gene_id = col_character(),
        gene_name = col_character(),
        baseMean = col_double(),
        log2FoldChange = col_double(),
        lfcSE = col_double(),
        stat = col_double(),
        pvalue = col_double(),
        padj = col_double()
      )
    )
    result_tables[[result_name]] <- de_table
  }
  de_results[[deseq2_setting[["name"]]]] <- result_tables
}

message("Gene enrichment online...")

# Individual
analysis_names <- purrr::map(settings[["gene_enrichment"]], "name")
gsea_outdir <- file.path(settings[["output_directory"]], "gene_set_enrichment")
dir.create(gsea_outdir, showWarnings = FALSE, recursive = TRUE)
gsea_enrichments <- purrr::map(settings[["gene_enrichment"]], function(settings_gene_enrichment) {
  analysis_name <- settings_gene_enrichment[["name"]]
  dir.create(file.path(gsea_outdir, analysis_name), showWarnings = FALSE, recursive = TRUE)
  results_names <- purrr::map(settings_gene_enrichment[["results"]], "name")
  out <- purrr::map(settings_gene_enrichment[["results"]], function(result) {
    result_name <- result[["name"]]
    message(sprintf(" - %s-%s", analysis_name, result_name))
    de_table <- de_results[[analysis_name]][[result_name]]
    if (is.null(de_table)) {
      print(names(de_results))
      print(names(de_results[[analysis_name]]))
      stop("de_table cant be null")
    }
    geneList <- de_table[["log2FoldChange"]]
    names(geneList) <- remove_ensembl_version(de_table[["gene_id"]])
    geneList <- sort(geneList, decreasing = TRUE)
    orgdb_obj <- AnnotationDbi::get(org_db)


    gsea_results <- c("BP", "MF", "CC") %>%
      purrr::set_names() %>%
      purrr::map(function(ont) {
        message("    - calling gseGO for ", length(geneList), " genes, on ", ont, " ontology")
        res <- gseGO_memoise(
          cache_dirname = settings[["cache"]],
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
        message("    - setting readable")
        res <- clusterProfiler::setReadable(res, org_db, 'ENSEMBL')
        message("    - term similarity")
        res <- tryCatch(
          {
            enrichplot::pairwise_termsim(res)
          },
          error = function(e) res
        )

        message("    - plotting")
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
              gsea_outdir,
              sprintf("%s-%s-%s_dotplot.png", analysis_name, result_name, ont)
            ),
            width = 20,
            height = 25,
            units="cm"
          )
        }
        write.csv(
          as.data.frame(res),
          file = file.path(gsea_outdir, analysis_name, sprintf("%s-%s.csv", result_name, ont)),
          row.names = FALSE
        )
        res
      })
    saveRDS(gsea_results, file.path(gsea_outdir, analysis_name, sprintf("%s-gseaGO.rds", result_name)))
    gsea_results
  })
  names(out) <- results_names
  out
})
names(gsea_enrichments) <- analysis_names
message("DONE")

