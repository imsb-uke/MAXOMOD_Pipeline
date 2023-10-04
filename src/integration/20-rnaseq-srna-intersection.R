#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/remove_ensembl_version.R"))

library(dplyr)

settings <- get_settings(omic = "integration", stage = "rnaseq_srna_intersection", dataset_name = NULL) # dataset_name="c9"

message("Reading srna de results with targets...")
de_results_thresholds_with_targets <- readRDS(file.path(settings[["srna_deg_digestion_dir"]], "de_results_with_targets.rds"))

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

message("Preparing miRNA & rnaseq intersection...")

message("Reading rnaseq results...")
rnaseq_de_result_thres <- readRDS(
  file.path(settings[["rnaseq_de_digested_dir"]], "de_tables", "de_result_thresholds.rds")
)


message("Reading ENSG to symbol mapping...")
ensg_symbol <- readRDS(file.path(settings[["rnaseq_prefiltpca_dir"]], "ENSG_to_symbol.rds"))
ensg_symbol$gene_id <- remove_ensembl_version(ensg_symbol$gene_id)
ensembl_to_symbol <- tibble::deframe(ensg_symbol)

dir.create(
  file.path(settings[["output_directory"]], "rnaseq_srna_intersection"),
  showWarnings = FALSE,
  recursive = TRUE
)

rnaseq_srna_intersect <- list()

for (analysis_name in intersect(names(rnaseq_de_result_thres), names(de_results_thresholds_with_targets$mature))) {
  rnaseq_srna_intersect[[analysis_name]] <- list()
  rna_analysis <- rnaseq_de_result_thres[[analysis_name]]
  srna_analysis <- de_results_thresholds_with_targets$mature[[analysis_name]]
  for (result_name in intersect(names(rna_analysis), names(srna_analysis))) {
    rnaseq_srna_intersect[[analysis_name]][[result_name]] <- list()
    rna_result <- rna_analysis[[result_name]]
    srna_result <- srna_analysis[[result_name]]
    for (criterion in intersect(names(rna_result), names(srna_result))) {
      rna_table <- rna_result[[criterion]]
      srna_table <- srna_result[[criterion]]
      rna_table$gene_id <- remove_ensembl_version(rna_table$gene_id)
      rnaseq_up_gene_id <- rna_table$gene_id[rna_table$regulation == "upregulated"]
      rnaseq_down_gene_id <- rna_table$gene_id[rna_table$regulation == "downregulated"]
      out <- srna_table %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          rnaseq_upregulated = list(intersect(ensembl_gene_id,  rnaseq_up_gene_id)),
          rnaseq_downregulated = list(intersect(ensembl_gene_id,  rnaseq_down_gene_id)),
          num_rnaseq_upregulated = length(rnaseq_upregulated),
          num_rnaseq_downregulated = length(rnaseq_downregulated)
        ) %>%
        dplyr::ungroup()
      out <- dplyr::filter(out, num_rnaseq_upregulated > 0 | num_rnaseq_downregulated > 0)
      out$rnaseq_upregulated <- purrr::map_chr(out$rnaseq_upregulated, ~ paste0(ensembl_to_symbol[.], collapse= ", "))
      out$rnaseq_downregulated <- purrr::map_chr(out$rnaseq_downregulated, ~ paste0(ensembl_to_symbol[.], collapse= ", "))
      out_print <- dplyr::select(
        out,
        dplyr::one_of(c("mirnaID", "baseMean", "log2FoldChange", "padj",
                        "mirna_regulation" = "regulation", "num_targets", "num_rnaseq_upregulated",
                        "num_rnaseq_downregulated", "rnaseq_upregulated", "rnaseq_downregulated")))
      write.csv(
        x = out_print,
        file = file.path(settings[["output_directory"]], "rnaseq_srna_intersection", sprintf("%s-%s-%s.csv", analysis_name, result_name, criterion)),
        row.names = FALSE
      )
      saveRDS(
        out,
        file = file.path(settings[["output_directory"]], "rnaseq_srna_intersection", sprintf("%s-%s-%s.rds", analysis_name, result_name, criterion))
      )
      rnaseq_srna_intersect[[analysis_name]][[result_name]][[criterion]] <- out
    }
  }
}

saveRDS(rnaseq_srna_intersect, file = file.path(settings[["output_directory"]], "rnaseq_srna_intersect.rds"))

message("DONE")
