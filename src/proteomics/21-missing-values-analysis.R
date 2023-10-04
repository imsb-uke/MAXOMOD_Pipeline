#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("naniar")
  library("readr")
  library("tibble")
  library("tidyr")
})


models <- tibble::tribble(
  ~name, ~dir_name,
  "sod1", "SOD1-mouse",
  "tdp43", "TDP43-mouse",
  "c9", "C9orf72-mouse",
  "fus", "FUS-mouse",
  "human", "human-datasets"
)

nainfo_per_model <- tibble::deframe(models) %>%
  purrr::map(
    function(model_dirname) {
      proteomics_base_dir <- sprintf("datasets/consortium/%s/02_organized_data/proteomics", model_dirname)
      two_matrix_files <- list(
        "orig" = list(
          file = file.path(proteomics_base_dir, "organized/intensity.rds"),
          dim_samples = "cols"
        ),
        "filtered" = list(
          file = file.path(proteomics_base_dir, "prefiltering-pca/intensity_mat_filtered.rds"),
          dim_samples = "rows"
        )
      )
      purrr::map(two_matrix_files, function(file_and_dim) {
        mat <- readRDS(file_and_dim$file)
        mat[mat == 0] <- NA
        if (file_and_dim$dim_samples == "cols") {
          mat <- t(mat)
        }
        num_samples <- nrow(mat)
        num_na_per_protein <- colSums(is.na(mat))
        table_num_na_per_protein <- table(num_na_per_protein)
        list(
          intensity = mat,
          num_samples = num_samples,
          num_na_per_protein = num_na_per_protein,
          table_num_na_per_protein = table_num_na_per_protein
        )
      }
      )
    }
  )


missingness_distr_combined <- 
  purrr::map_dfr(nainfo_per_model, function(model) {
    purrr::map_dfr(model, function(x) {
      as.data.frame(x$table_num_na_per_protein) %>%
        mutate(num_na_per_protein = as.numeric(as.character(num_na_per_protein))) %>%
        arrange(num_na_per_protein) %>%
        mutate(CumFreq = cumsum(Freq),
               CumFrac = CumFreq/sum(Freq),
               frac_na_per_protein = num_na_per_protein/x$num_samples
        )
    },
    .id = "step"
    )
  },
  .id = "dataset"
  ) %>%
  tibble::as_tibble() %>%
  mutate(
    dataset = factor(dataset, levels = models$name),
    step = factor(step, levels = c("orig", "filtered"))
  ) %>%
  dplyr::select(dataset, step, num_na_per_protein, frac_na_per_protein, Freq, CumFreq, CumFrac)


gplt <- ggplot(missingness_distr_combined) +
  geom_col(aes(x = num_na_per_protein, y = Freq)) +
  labs(
    y = "<y> number of proteins", 
    x = "missing in <x> samples",
    title = "Missing values of proteins") +
  facet_grid(step~dataset, scales = "free")

dir.create("results/proteomics-missing-values", recursive = TRUE, showWarnings = FALSE)
ggsave(filename = "results/proteomics-missing-values/missing-value-dist.png", plot = gplt)
saveRDS(missingness_distr_combined, file = "results/proteomics-missing-values/missing-value-dist.rds")

gplt2 <- ggplot(missingness_distr_combined) +
  geom_line(aes(y = frac_na_per_protein, x = CumFrac, color = step)) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "This percentage of proteins...", 
    y = "...has this or less samples missing",
    title = "Missing values of proteins") +
  facet_grid(~dataset)
