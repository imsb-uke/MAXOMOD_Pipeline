#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("readr")
  library("tibble")
})

settings <- get_settings(omic = "phosphoproteomics", stage = "organize_samples", dataset_name = NULL) # dataset_name = "c9"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = cols()
)

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)


received_proteomics_data <- readr::read_tsv(settings[["maxquant_input"]], col_types = cols())


received_proteomics_data_simpleIDs <- received_proteomics_data %>%
    mutate(id_first = gsub(";.*", "", `Proteins`)) %>%
    tidyr::separate(id_first, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|") %>%
    dplyr::filter(Position!='Non numérique') %>%
    dplyr::mutate(PhosphositeName=paste0(Position,`Amino acid`,'_',UniProtName,Multiplicity))

    #dplyr::mutate(UID_UniProtName=paste(`Position`,UniProtName,sep='_')) %>%
    #distinct(UID_UniProtName, .keep_all = TRUE)

#protein_group_to_proteins <- received_proteomics_data %>%
#  dplyr::select(`Proteins`) %>%
#  mutate(id_first = gsub(";.*", "", `Proteins`)) %>%
#  tidyr::separate(id_first, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|") %>%
#  dplyr::mutate(AllIds = strsplit(`Proteins`, split = ";", fixed = TRUE)) %>%
#  dplyr::select(UniProtGroup = UniProtName, AllIds) %>%
#  tidyr::unnest(AllIds) %>%
#  tidyr::separate(AllIds, into = c("SwissProt_Tremble", "UniProtAccession", "UniProtName"), sep = "\\|") %>%
#  dplyr::filter(SwissProt_Tremble == "sp", !is.na(UniProtGroup)) %>%
#  dplyr::select(UniProtGroup, UniProtAccession, UniProtName) %>%
#  dplyr::group_by(UniProtGroup) %>%
#  dplyr::mutate(NumProteinsinGroup = n()) %>%
#  dplyr::ungroup()


intensity <- received_proteomics_data_simpleIDs %>%
  dplyr::filter(SwissProt_Tremble == "sp", UniProtName != "Biognosys_iRT") %>%
  dplyr::select(PhosphositeName, starts_with(settings[["intensity_column_prefix"]])) %>%
  dplyr::rename_with(~ gsub(settings[["intensity_name_pattern"]], settings[["intensity_name_replacement"]], .),
                     starts_with(settings[["intensity_column_prefix"]]) )


if (!is.null(settings[["rename_samples"]])) {
  sample_map <- readr::read_tsv(settings[["rename_samples"]], col_types = cols()) %>%
    dplyr::select(2,1) %>%
    tibble::deframe()
  if (length(sample_map) > 0) {
    print(colnames(intensity))
    intensity <- intensity %>%
      dplyr::rename(!!!sample_map)
    print(colnames(intensity))
  }
}



intensity <- intensity %>%
  dplyr::select(PhosphositeName, one_of(cohort$SampleID))

#NA in phosphoproteomics data needs to be replaced by an intensity of 0
intensity <- intensity %>%
  replace(is.na(.), 0)

#i_mat <- intensity %>%
#  tibble::column_to_rownames(UID_UniProtName) %>%
#  as.matrix()
#protein_group_to_proteins <- left_join(
#  protein_group_to_proteins,
#  tibble::enframe(rowSums(i_mat == 0), name = "UniProtGroup", value = "SamplesMissing"),
#  by = "UniProtGroup"
#)



#if (!all(c("UID_UniProtName", cohort$SampleID) == colnames(intensity))) {
#  print("Cohort sample IDs:")
#  print(cohort$SampleID)
#  print("Intensity matrix ids")
#  print(colnames(intensity)[-1])
#  stop("Mismatch proteomics samples and cohort")
#}


protein_annotations <- received_proteomics_data_simpleIDs %>%
    dplyr::filter(SwissProt_Tremble == "sp", UniProtName != "Biognosys_iRT") %>%
    dplyr::select(PhosphositeName,UniProtName, UniProtAccession,
                  `Amino acid`,`Charge`,`Reverse`,`Potential contaminant`,`Protein group IDs`,`Fasta headers`,`Multiplicity`,
                  `Localization prob`,PEP,Score,`Delta score`,`Score for localization`,`Mass error [ppm]`,Intensity,Position,Proteins,`Positions within proteins`,`Leading proteins`,Protein,`Sequence window`,`Phospho (STY) Probabilities`,id,`Unique identifier`,
                  `Number of Phospho (STY)`, SwissProt_Tremble)


write.csv(x = intensity, file = file.path(settings[["output_directory"]], "intensity.csv"), row.names = FALSE)
write.csv(x = protein_annotations, file = file.path(settings[["output_directory"]], "protein_annotations.csv"), row.names = FALSE)
#write.csv(x = protein_group_to_proteins, file = file.path(settings[["output_directory"]], "protein_group_to_proteins.csv"), row.names = FALSE)

saveRDS(intensity, file = file.path(settings[["output_directory"]], "intensity.rds"))
saveRDS(protein_annotations, file = file.path(settings[["output_directory"]], "protein_annotations.rds"))
#saveRDS(protein_group_to_proteins, file = file.path(settings[["output_directory"]], "protein_group_to_proteins.rds"))
