#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("readr")
  library("tibble")
})

settings <- get_settings(omic = "phosphoproteomics", stage = "filter_STY", dataset_name = NULL) # dataset_name = "c9"
loc_prob_thresh <- settings[["filter_loc_prob"]] 

raw_STY <- readr::read_tsv(settings[["input"]], col_types = cols())

filtered_STY <- raw_STY %>%
    dplyr::filter(is.na(`Potential contaminant`), #removing possible contaminants -> '+' if possible contaminant, otherwise NA
                    is.na(Reverse), #removing reverse proteins -> '+' if reverse, otherwise NA
                    `Localization prob` >= loc_prob_thresh)

message(paste("Reduced Data from",nrow(raw_STY),"to",nrow(filtered_STY),"rows...",sep=' '))
message(paste("Now without contaminants, reverse peptides and peptides with localization probability of below",loc_prob_thresh,"...",sep=' '))

message("Reducing col number...")
reduced_STY <- filtered_STY %>% 
    dplyr::select(matches(settings[["intensity_name_pattern"]]),
                  any_of(settings[["keep_cols"]]))

message(paste("Saving filtered and reduced file to",settings[["output"]],"...",sep=' '))
write_delim(reduced_STY,settings[["output"]],delim='\t')

message("DONE")
