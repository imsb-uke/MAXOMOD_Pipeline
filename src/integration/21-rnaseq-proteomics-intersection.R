#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  requireNamespace("gtools")
})

source(here::here("src/shared/R-helpers/remove_ensembl_version.R"))
source(here::here("src/shared/R-helpers/get_settings.R"))


settings <- get_settings(omic = "integration", stage = "rnaseq_proteomics_intersection", dataset_name = NULL) # dataset_name="sod1"

# settings$output_directory <- "results/integration/rnaseq_proteomics/SOD1-mouse/"
# settings$uniprot_to_ensembl <- "database/uniprot/uniprot_MOUSE_10090_idmapping_selected.tab.gz"
# settings$rnaseq_de_results <- "datasets/consortium/SOD1-mouse/02_organized_data/rnaseq/deg-digestion/de_tables/covariate_sex-mut_vs_ctrl-strict_justpval.csv"
# settings$protein_group_to_proteins <- "datasets/consortium/SOD1-mouse/02_organized_data/proteomics/organized/protein_group_to_proteins.rds"
# settings$proteomics_de_results <- "datasets/consortium/SOD1-mouse/02_organized_data/proteomics/deg/covariate_sex/mut_vs_ctrl.rds"

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

uniprot_to_ensembl <- readr::read_tsv(settings[["uniprot_to_ensembl"]], col_names = FALSE)
uniprot_to_ensembl <- uniprot_to_ensembl[, c("X2", "X19")]
colnames(uniprot_to_ensembl) <- c("UniProtName", "ENSEMBL")


rnaseq_de_results <- readr::read_csv(settings[["rnaseq_de_results"]], col_types = cols())
protein_group_to_proteins <- readRDS(settings[["protein_group_to_proteins"]])
prot_de_results <- readRDS(settings[["proteomics_de_results"]])

num_rnaseq_deg <- nrow(rnaseq_de_results)

rnaseq_deg_with_candidate_proteins <- rnaseq_de_results %>%
  dplyr::mutate(ENSEMBL = remove_ensembl_version(gene_id)) %>%
  dplyr::left_join(uniprot_to_ensembl, by = "ENSEMBL") %>%
  dplyr::mutate(MeasuredInProteomics = UniProtName %in% protein_group_to_proteins$UniProtName) %>%
  dplyr::filter(MeasuredInProteomics) %>%
  dplyr::left_join(
    protein_group_to_proteins %>% dplyr::select(UniProtName, UniProtGroup, NumProteinsinGroup, SamplesMissing),
    by = "UniProtName"
  ) %>%
  dplyr::left_join(
    prot_de_results %>% dplyr::rename(UniProtGroup = UniProtName),
    by = "UniProtGroup",
    suffix = c("_rna", "_prot")
  ) %>%
  dplyr::filter(!is.na(baseMean_prot))

num_rnaseq_deg_with_detected_protein <- nrow(rnaseq_deg_with_candidate_proteins)

num_rnaseq_deg_with_signif_protein <- nrow(dplyr::filter(rnaseq_deg_with_candidate_proteins, padj_prot < 0.05))

info_msg <- sprintf(
  "There are %d differentially expressed genes, %d have a detected protein, and %d have a differentially abundant protein",
  num_rnaseq_deg,
  num_rnaseq_deg_with_detected_protein,
  num_rnaseq_deg_with_signif_protein
)

message(info_msg)
writeLines(info_msg, file.path(settings[["output_directory"]], "info.txt"))

rnaseq_deg_with_candidate_proteins_table <- rnaseq_deg_with_candidate_proteins %>% 
  transmute(
    gene = gene_name,
    baseMean_rna,
    log2FC_rna = log2FoldChange_rna,
    pv_rna = gtools::stars.pval(pvalue_rna),
    padj_rna = gtools::stars.pval(padj_rna),
    prot_name = UniProtName,
    nprot_group = NumProteinsinGroup, 
    SampMiss = SamplesMissing,
    baseMean_prot,
    log2FC_prot = log2FoldChange_prot,
    pv_prot = gtools::stars.pval(pvalue_prot), 
    padj_prot = gtools::stars.pval(padj_prot)
  )

write.csv(
  x = rnaseq_deg_with_candidate_proteins_table,
  file = file.path(settings[["output_directory"]], "rnaseq_and_proteomics.csv"),
  row.names = FALSE
)

#rnaseq_deg_with_candidate_proteins_table %>% View(title = "RNASeq & Proteomics")

gplt <- ggplot(rnaseq_deg_with_candidate_proteins) + 
  geom_point(aes(x = log2FoldChange_rna, y = log2FoldChange_prot, color = padj_prot, label = gene_name)) +
  scale_color_gradient2(limits = c(0, 0.1), low = "blue", mid = "black", high = "black", midpoint =0.1) +
  labs(
    x = "log2FC (rnaseq)",
    y = "log2FC (proteomics)",
    color = "p.adj (prot)",
    title = "RNASeq and proteomics log2 fold changes"
  )

ggsave(
  filename = file.path(settings[["output_directory"]], "rnaseq_and_proteomics.png"),
  plot = gplt,
  width = 12,
  height = 12,
  unit = "in"
)

#ggplotly(gplt)
message("DONE")

