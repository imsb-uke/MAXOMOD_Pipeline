#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("readr")
  library("tibble")
  library("ggplot2")
  requireNamespace("forcats")
  library("preprocessCore")
  library("protti") # framework used for QC plots
  library("data.table")
})


settings <- get_settings(omic = "proteomics", stage = "additional_QC", dataset_name = NULL) #dataset_name = "sod1"


cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = NULL
) %>% dplyr::filter(!SampleID %in% settings[["exclude_samples"]])


print(cohort)

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

#protein_intensity <- readRDS(file.path(settings[["input_directory"]], "intensity.rds"))
#protein_annotations <- readRDS(file.path(settings[["input_directory"]], "protein_annotations.rds"))

protein_data <-readr::read_csv(settings[["input_data"]])

print(str(protein_data))
print(head(protein_data))

protein_longer <- protein_data %>%
  pivot_longer(
    cols = !settings[["pivot_longer_exclude_col"]],
    names_to = settings[["pivot_longer_names_to"]],
    values_to = settings[["pivot_longer_values_to"]]
  ) %>%
  dplyr::rename(!!settings[["rename_index_new"]] := !!settings[["rename_index_old"]]) #rename index col

print(str(protein_longer))
print(head(protein_longer))


#transforms back log2 values for e.g. Coefficients of variation QC plats
if (settings[["isLog2"]]){
    protein_longer <- protein_longer %>%
        mutate(raw_intensity = as.integer(2^intensity))
}

print(str(protein_longer))
print(head(protein_longer))

# Combine your long data frame with the sample annotation
protein_longer_annotated <- protein_longer %>% 
  left_join(y = cohort, by = "SampleID")

print(str(protein_longer_annotated))
print(head(protein_longer_annotated))

qc_cvs(data = protein_longer_annotated,
       grouping = proteinID,
       condition = Condition,
       intensity = !!as.symbol(ifelse(settings[["isLog2"]],"raw_intensity",settings[["pivot_longer_values_to"]])), 
       plot = TRUE, 
       plot_style = "violin")
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_cvs.pdf")))

protein_longer_annotated$intensity <- replace(protein_longer_annotated$intensity, protein_longer_annotated$intensity == 0, NA)
    
qc_ids(data = protein_longer_annotated,
       sample = SampleID,
       grouping = proteinID,
       intensity = intensity,
       condition = Condition, 
       title = "Protein identifications per sample",
       plot = TRUE)
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_protein_identifications_per_sample.pdf")),width=settings[["qc_ids_width"]])
              
qc_protein_idenfications <- qc_ids(data = protein_longer_annotated,
       sample = SampleID,
       grouping = proteinID,
       intensity = intensity,
       condition = Condition, 
       title = "Protein identifications per sample",
       plot = FALSE)
write_csv(qc_protein_idenfications,
          file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_protein_identifications_per_sample.csv")))

qc_median_intensities(data = protein_longer_annotated,
                      sample = SampleID,
                      grouping = proteinID,
                      intensity = intensity)
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_median_intensities.pdf")),width=settings[["qc_ids_width"]])

qc_data_completeness(data = protein_longer_annotated, 
                     sample = SampleID, 
                     grouping = proteinID, 
                     intensity = intensity, 
                     plot = TRUE)
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_data_completeness.pdf")),width=settings[["qc_dc_width"]])

if(settings[["isLog2"]]){
    qc_intensity_distribution(data = protein_longer_annotated,
                          sample = SampleID,
                          grouping = proteinID,
                          intensity_log2 = intensity,
                          plot_style = "boxplot")
    ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_intensity_distribution_boxplot.pdf")),width=settings[["qc_ids_width"]])
    
    qc_intensity_distribution(data = protein_longer_annotated,
                          grouping = proteinID,
                          intensity_log2 = intensity,
                          plot_style = "histogram")
    ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_intensity_distribution_histogram.pdf")),width=settings[["qc_ids_width"]])

sample_correlation_plot <- qc_sample_correlation(data = protein_longer_annotated,
                      sample = SampleID,
                      grouping = proteinID,
                      intensity_log2 = intensity, 
                      condition = Condition, 
                      interactive = FALSE)
 ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_sample_correlation.pdf")),plot=sample_correlation_plot,width=settings[["qc_ids_width"]],height=settings[["qc_ids_width"]])
}

qc_pca(
  data = protein_longer_annotated,
  sample = SampleID,
  grouping = proteinID,
  intensity = intensity,
  condition = Condition,
  digestion = NULL,
  plot_style = "scree"
)
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_pca_scree.pdf")),width=settings[["qc_pca_width"]])

qc_pca(
  data = protein_longer_annotated,
  sample = SampleID,
  grouping = proteinID,
  intensity = intensity,
  condition = Condition,
  components = c("PC1", "PC2"), 
  plot_style = "pca"
)
ggsave(file.path(settings[["output_directory"]],paste0(unique(protein_longer_annotated$CellLine),"_qc_pca.pdf")),width=settings[["qc_ids_width"]])