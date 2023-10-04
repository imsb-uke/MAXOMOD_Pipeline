#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/memoise_decorator.R"))
source(here::here("src/shared/R-helpers/enrichGO_memoise.R"))
source(here::here("src/shared/R-helpers/remove_ensembl_version.R"))

suppressPackageStartupMessages({
  library("argparse")
  library("cowplot")
  library("dplyr")
  library("DESeq2")
  library("forcats")
  library("ggplot2")
  library("matrixStats")
  library("pheatmap")
  library("purrr")
  library("RColorBrewer")
  library("readr")
  library("tibble")
  library("tidyr")
  library("yaml")
  library("enrichplot")
})

settings <- get_settings(omic="rnaseq", stage="gene_gene_correlation", dataset_name=NULL) #dataset_name = "c9")

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

if (settings[["organism"]] == "mmusculus") {
  org_db <- "org.Mm.eg.db"
} else if (settings[["organism"]] == "hsapiens") {
  org_db <- "org.Hs.eg.db"
} else {
  stop("Not implemented error")
}


# Gene gene correlation ------------------------
vsd <- readRDS(file = file.path(settings[["prefiltpca_dir"]], "counts_norm_vst.rds"))
rnaseq_norm_mat <- t(assay(vsd))

cluster_genes_into_modules_nomem <- function(rnaseq_norm_mat, power = 6) {
  connectivity <- WGCNA::softConnectivity(datExpr = rnaseq_norm_mat, power = power)

  # Plot a histogram of k and a scale free topology plot
  #hist(connectivity)
  #WGCNA::scaleFreePlot(connectivity, main="Check scale free topology\n")
  ADJ1 <- abs(cor(rnaseq_norm_mat,use="pairwise.complete.obs"))^power
  dissTOM <- WGCNA::TOMdist(ADJ1)

  hierTOM <- hclust(as.dist(dissTOM),method="average")
  labelsDynamicTOM <- dynamicTreeCut::cutreeDynamic(hierTOM, distM=dissTOM, method="hybrid")
  #colorDynamicTOM <- WGCNA::labels2colors(labelsDynamicTOM)
  names(labelsDynamicTOM) <- colnames(rnaseq_norm_mat)
  gene_id_to_module <- tibble::enframe(labelsDynamicTOM, name = "gene_id", value = "module")
  gene_id_to_module
}

cluster_genes_into_modules <- memoise_decorator(
  cache_dirname = settings[["cache"]],
  cache_prefix = "_cluster_genes_into_modules",
  fun = cluster_genes_into_modules_nomem
)

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = cols()
)

power <- settings[["power"]]
gene_id_to_module <- list(
  all = cluster_genes_into_modules(rnaseq_norm_mat, power = power),
  females = cluster_genes_into_modules(rnaseq_norm_mat[cohort$SampleID[cohort$Sex == "female"],], power = power),
  males = cluster_genes_into_modules(rnaseq_norm_mat[cohort$SampleID[cohort$Sex == "male"],], power = power)
)

saveRDS(gene_id_to_module, file.path(settings[["output_directory"]], "gene_id_to_module.rds"))

named_modules <- sort(unique(gene_id_to_module$all$module))
names(named_modules) <- sprintf("Module_%02d", named_modules)

module_enrichment <- purrr::map(
  named_modules,
  function(mod) {
    purrr::map(
      c("BP" = "BP", "CC" = "CC", "MF" = "MF", "ALL" = "ALL"),
      function(ont) {
        enrichGO_memoise(
          cache_dirname = settings[["cache"]],
          gene = gene_id_to_module$all %>%
            dplyr::filter(module == !!mod) %>%
            dplyr::pull("gene_id") %>%
            remove_ensembl_version(),
          OrgDb = org_db,
          universe = remove_ensembl_version(gene_id_to_module$all$gene_id),
          keyType = "ENSEMBL",
          ont = ont,
          minGSSize = 10,
          maxGSSize = 500,
          readable = TRUE
        )
    }
    )
  }
)

saveRDS(module_enrichment, file.path(settings[["output_directory"]], "module_enrichment.rds"))

module_plot_outdir <- file.path(settings[["output_directory"]], "module_plot_strict")
dir.create(module_plot_outdir, recursive = TRUE, showWarnings = FALSE)

message("Simplifying module enrichment results...")
module_enrichment_simplified <- purrr::map(names(module_enrichment), function(mod_list_name) {
  mod_list <- module_enrichment[[mod_list_name]]
  purrr::map(c("BP" = "BP", "CC" = "CC", "MF" = "MF"), function(ont) {
    message(paste0(mod_list_name, "-", ont))
    mod <- mod_list[[ont]]
    if (nrow(as.data.frame(mod)) == 0) {
      return(NULL)
    }
    message("...simplify")
    xx <- clusterProfiler::simplify(mod)
    message("...pairwise_termsim")
    xx <- enrichplot::pairwise_termsim(xx)
    write.csv(as.data.frame(xx),
              file = file.path(module_plot_outdir, paste0(mod_list_name, "_", ont, "_parent_terms.csv")))
    message("...plots")
    tryCatch({
      plt <- enrichplot::emapplot(xx, layout = "kk")
      ggplot2::ggsave(
        filename = file.path(module_plot_outdir, paste0(mod_list_name, "_", ont, "_parent_terms.pdf")),
        plot = plt,
        width = 15,
        height = 5,
        units = "in"
      )
      ggplot2::ggsave(
        filename = file.path(module_plot_outdir, paste0(mod_list_name, "_", ont, "_parent_terms.png")),
        plot = plt,
        width = 15,
        height = 5,
        units = "in"
      )
    },
    error=function(cond) {
    })
    xx
  })
})

saveRDS(module_enrichment_simplified, file.path(settings[["output_directory"]], "module_enrichment_simplified.rds"))

message("Creating _dotplot for each module and ontology")
purrr::map(names(module_enrichment), function(mod_list_name) {
  mod_list <- module_enrichment[[mod_list_name]]
  purrr::map(c("BP" = "BP", "CC" = "CC", "MF" = "MF", "ALL" = "ALL"), function(ont) {
    message(paste0(mod_list_name, "-", ont))
    mod <- mod_list[[ont]]
    if (nrow(as.data.frame(mod)) == 0) {
      return(NULL)
    }
    plt <- dotplot(mod, showCategory=20) +
      ggplot2::ggtitle(paste0(mod_list_name, "-", ont))

    ggplot2::ggsave(
      filename = file.path(module_plot_outdir, paste0(mod_list_name, "_" , ont, "_dotplot.pdf")),
      plot = plt,
      width = 15,
      height = 5,
      units = "in"
    )

    ggplot2::ggsave(
      filename = file.path(module_plot_outdir, paste0(mod_list_name, "_" , ont, "_dotplot.png")),
      plot = plt,
      width = 15,
      height = 5,
      units = "in"
    )
  })
})



fold_change_as_vec <- function(df, padj_thres = 0.05, log2fc_cap = 2) {
  df %>%
    dplyr::filter(padj < !!padj_thres) %>%
    dplyr::transmute(
      gene_id = remove_ensembl_version(gene_id),
      log2FoldChange=dplyr::if_else(
        abs(log2FoldChange) > log2fc_cap,
        sign(log2FoldChange)*log2fc_cap,
        log2FoldChange
      )
    ) %>%
    tibble::deframe()
}


message("Reading DEG outputs")

mut_vs_ctrl_fem <- readr::read_csv(
  file = file.path(settings[["deg_dir"]], "only_females", "mut_vs_ctrl.csv"),
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


mut_vs_ctrl_male <- readr::read_csv(
  file = file.path(settings[["deg_dir"]], "only_males", "mut_vs_ctrl.csv"),
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


module_plot_outdir <- file.path(settings[["output_directory"]], "module_plot_strict")
dir.create(module_plot_outdir, recursive = TRUE, showWarnings = FALSE)

module_netplot <- purrr::map(
  names(module_enrichment),
  function(mod_name, padj_thres) {
    message(mod_name)
    module <- module_enrichment[[mod_name]][["ALL"]]
    if (nrow(as.data.frame(module)) == 0) {
      return(NULL)
    }
    p1 <- cnetplot(
      module,
      foldChange = fold_change_as_vec(mut_vs_ctrl_fem, padj_thres = padj_thres),
      cex_label_category = 0.5,
      cex_label_gene = 0.5
    ) + ggplot2::ggtitle(label = mod_name, subtitle = "Mut vs ctrl (females)")

    p2 <- cnetplot(
      module,
      foldChange = fold_change_as_vec(mut_vs_ctrl_male, padj_thres = padj_thres),
      cex_label_category = 0.5,
      cex_label_gene = 0.5
    )+ ggplot2::ggtitle(label = mod_name, subtitle = "Mut vs ctrl (males)")

    out <- cowplot::plot_grid(p1, p2, ncol = 2)
    saveRDS(out, file = file.path(module_plot_outdir, paste0(mod_name, "_netplot.rds")))
    filename <- file.path(module_plot_outdir, paste0(mod_name, "_netplot.pdf"))
    cowplot::save_plot(
      filename = filename,
      plot = out,
      base_height = 15
    )
    filename <- file.path(module_plot_outdir, paste0(mod_name, "_netplot.png"))
    cowplot::save_plot(
      filename = filename,
      plot = out,
      base_height = 15
    )
    out
  },
  padj_thres = 0.05
)


#
#
#
# # strict
#
mut_vs_ctrl_male %>%
  dplyr::left_join(gene_id_to_module$all, by = "gene_id") %>%
  dplyr::mutate(
    is_significant = !is.na(padj) & padj < 0.05,
    is_upreg = is_significant & log2FoldChange > 0,
    is_downreg = is_significant & log2FoldChange < 0
  ) %>%
  group_by(module) %>%
  dplyr::summarize(
    total = n(),
    num_up = sum(is_upreg),
    num_down = sum(is_downreg),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc((num_up + num_down)/total))
#
#
compute_fisher_test_module_enrichment <- function(de_result, gene_id_to_module) {
  total_signif_genes <- sum(!is.na(de_result$padj) & de_result$padj < 0.05)
  total_genes <- nrow(de_result)
  out <- de_result %>%
    dplyr::left_join(gene_id_to_module, by = "gene_id") %>%
    dplyr::mutate(
      is_significant = !is.na(padj) & padj < 0.05,
    ) %>%
    dplyr::count(module, is_significant) %>%
    tidyr::spread(is_significant, n, fill=0) %>%
    dplyr::transmute(
      module = module,
      n_not_signif = `FALSE`,
      n_signif = `TRUE`,
      frac_signif = 100*`TRUE`/(`FALSE` + `TRUE`)
    ) %>%
    dplyr::arrange(desc(frac_signif))

  out$odds_ratio <- NA_real_
  out$odds_ratio_pval <- NA_real_
  for (i in seq_len(nrow(out))) {
    xx <- fisher.test(
      matrix(
        data = c(out$n_signif[i], total_signif_genes,
                 out$n_signif[i] + out$n_not_signif[i], total_genes),
        nrow = 2,
        ncol = 2,
        byrow = TRUE
      )
    )
    out$odds_ratio[i] <- xx$estimate
    out$odds_ratio_pval[i] <- xx$p.value
  }
  out$odds_ratio_padj <- p.adjust(out$odds_ratio_pval, method = "BH")
  out$total_signif_genes <- total_signif_genes
  out$total_tested_genes <- total_genes
  out$total_proportion <- total_signif_genes/total_genes
  out
}

mut_vs_ctrl_male_with_modules <- compute_fisher_test_module_enrichment(de_result=mut_vs_ctrl_male, gene_id_to_module=gene_id_to_module$all)
mut_vs_ctrl_female_with_modules <- compute_fisher_test_module_enrichment(de_result=mut_vs_ctrl_fem, gene_id_to_module=gene_id_to_module$all)

mut_vs_ctrl_with_modules_all <- bind_rows(
  male = mut_vs_ctrl_male_with_modules ,
  female = mut_vs_ctrl_female_with_modules,
  .id = "sex"
)

saveRDS(
  mut_vs_ctrl_with_modules_all,
  file = file.path(module_plot_outdir, "mut_vs_ctrl_with_modules_all.rds")
)

mut_vs_ctrl_with_modules_high <- mut_vs_ctrl_with_modules_all %>%
  dplyr::filter(odds_ratio > 1, odds_ratio_padj < 0.05)

saveRDS(
  mut_vs_ctrl_with_modules_high,
  file = file.path(module_plot_outdir, "mut_vs_ctrl_with_modules_high.rds")
)



# Module enrichment, now only the differentially expressed within the module:
compute_go_enrich_on_diffexpr_per_module <- function(de_result, gene_id_to_module, named_modules) {
  de_results_with_modules <- de_result %>%
    dplyr::left_join(gene_id_to_module, by = "gene_id") %>%
    dplyr::mutate(
      is_significant = !is.na(padj) & padj < 0.05,
    )
  module_enrichment_difexpr <- purrr::map(
    named_modules,
    function(mod) {
      purrr::map(
        c("BP" = "BP", "CC" = "CC", "MF" = "MF", "ALL" = "ALL"),
        function(ont) {
          enrichGO_memoise(
            cache_dirname = settings[["cache"]],
            gene = de_results_with_modules %>%
              dplyr::filter(module == !!mod, is_significant == TRUE) %>%
              dplyr::pull("gene_id") %>%
              remove_ensembl_version(),
            OrgDb = org_db,
            universe = gene_id_to_module %>%
              dplyr::filter(module == !!mod) %>%
              dplyr::pull("gene_id") %>%
              remove_ensembl_version(),
            keyType = "ENSEMBL",
            ont = ont,
            minGSSize = 10,
            maxGSSize = 500,
            readable = TRUE
          )
        }
      )
    }
  )
  module_enrichment_difexpr
}

module_enrichment_diffexpr_male <- compute_go_enrich_on_diffexpr_per_module(mut_vs_ctrl_male, gene_id_to_module$all, named_modules)
module_enrichment_diffexpr_female <- compute_go_enrich_on_diffexpr_per_module(mut_vs_ctrl_fem, gene_id_to_module$all, named_modules)
saveRDS(module_enrichment, file.path(settings[["output_directory"]], "module_enrichment_diffexpr_male.rds"))
saveRDS(module_enrichment, file.path(settings[["output_directory"]], "module_enrichment_diffexpr_female.rds"))

