#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/create_volcano_plot.R"))

suppressPackageStartupMessages({
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
  library("EnhancedVolcano")
})

settings <- get_settings(omic = "srna", stage = "deg_digestion", dataset_name = NULL) # dataset_name="c9"

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

message("Reading DEG outputs")

get_de_results <- function(deseq2_settings, deg_dir) {
  de_results <- list()
  for (deseq2_setting in deseq2_settings) {
    deseq2_result_dir <- file.path(deg_dir, deseq2_setting[["name"]])
    result_tables <- list()
    for (result in deseq2_setting[["results"]]) {
      result_name <- result[["name"]]
      csv_file <- file.path(deseq2_result_dir, sprintf("%s.csv", result_name))
      de_table <- readr::read_csv(
        file = csv_file,
        col_types = cols(
          mirnaID = col_character(),
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
  de_results
}

de_results <- list()
for (mature_or_hairpin in c("mature", "hairpin")) {
  de_results[[mature_or_hairpin]] <- get_de_results(
    deseq2_settings = settings[["deseq2"]],
    deg_dir = file.path(settings[["deg_dir"]], mature_or_hairpin)
  )
}

volcano_single_plots_for_de_results <- function(de_results, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  volcano_plots <- purrr::imap(de_results, function(results, deseq2_setting_name) {
    volcano_results_dir <- file.path(outdir, deseq2_setting_name)
    dir.create(volcano_results_dir, showWarnings = FALSE, recursive = TRUE)
    purrr::imap(results, function(result, result_name) {
      plt <- create_volcano_plot(
        df = result,
        padj_thres = 0.05,
        up_log2fc_thres = log2(1.5),
        down_log2fc_thres = -log2(1.5),
        padj_thres_loose = 0.1,
        cap_minuslog10pval = 10,
        label_featurename = "mirnaID"
      )
      ggsave(filename = file.path(volcano_results_dir, sprintf("%s.pdf", result_name)),
             plot = plt, height = unit(4, "in"), width = unit(6, "in"))
        
      #create additional enhancedVolcano plots
      for (pval_thresh in c(0.05,0.1)) {
      enhancedplt<-EnhancedVolcano(result,
                    lab =result$mirnaID,
                    x = 'log2FoldChange',
                    y = 'padj',#'pvalue',#
                    pCutoff = pval_thresh,
                    FCcutoff = log2(1.5),
                    ylim = c(0, 
                             ifelse(-log10(pval_thresh)<(max(-log10(result$padj))+sd(-log10(result$padj))),
                                max(-log10(result$padj))+sd(-log10(result$padj)),
                                -log10(pval_thresh))),
                    xlim = c(min(result$log2FoldChange)-sd(result$log2FoldChange),
                             max(result$log2FoldChange)+sd(result$log2FoldChange)),
                    pointSize = 1,
                    title =sprintf('Volcano plot with a pvalue threshold of %s', pval_thresh),
                    subtitle ='NS = non significant',
                   )
      ggsave(filename = file.path(volcano_results_dir, paste0(result_name,"_enhanced",pval_thresh,".pdf")),
           plot = enhancedplt)#, height = unit(4, "in"), width = unit(6, "in"))
      }
      plt
    })
  })
  volcano_plots
}

volcano_merged_plots <- function(volcano_plot_objs, volcano_plot_settings, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (volcano_merged_plot in volcano_plot_settings) {
    name <- volcano_merged_plot[["name"]]
    subplots <- volcano_merged_plot[["subplots"]]
    plot_list <- list()
    for (subplot in subplots) {
      analysis <- subplot[["name"]][[1]]
      result <- subplot[["name"]][[2]]
      subplot_title <- subplot[["title"]]
      plt <- volcano_plot_objs[[analysis]][[result]]
      plt <- plt + labs(title = subplot_title)
      plot_list <- append(plot_list, list(plt))
    }
    plt <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
    ggsave(filename = file.path(outdir, sprintf("%s.pdf", name)),
           plot = plt, height = unit(6, "in"), width = unit(8, "in"))
    plt
  }
}

message("Creating volcano plots")
volcano_single_mature <- volcano_single_plots_for_de_results(
  de_results = de_results$mature,
  outdir = file.path(settings[["output_directory"]], "volcano_single", "mature")
)
volcano_merged_plots(
  volcano_plot_objs = volcano_single_mature,
  volcano_plot_settings = settings[["volcano_plots"]],
  outdir = file.path(settings[["output_directory"]], "volcano_merged", "mature")
)

volcano_single_hairpin <- volcano_single_plots_for_de_results(
  de_results = de_results$hairpin,
  outdir = file.path(settings[["output_directory"]], "volcano_single", "hairpin")
)
volcano_merged_plots(
  volcano_plot_objs = volcano_single_mature,
  volcano_plot_settings = settings[["volcano_plots"]],
  outdir = file.path(settings[["output_directory"]], "volcano_merged", "hairpin")
)

#Merge male and female volcano plot into one
for (srna_form in c("mature","hairpin")){
    for (pval_type in c('padj','pvalue')){
        for (pval_thresh in c(0.05,0.1)) { 
            results <- list()
            for (setting in settings[["merged_volcano"]]){
                analysis <- setting[["analysis"]]
                title <- setting[["title"]]
                result_name <- setting[["result_name"]]
                color <- setting[["color"]]
                result <- de_results[[srna_form]][[analysis]][[result_name]] %>%
                    mutate(Legend=if_else((between(log2FoldChange, -log2(1.5), log2(1.5))|(get(pval_type) >= pval_thresh)),"not significant",title),
                          Legend_color=if_else((between(log2FoldChange, -log2(1.5), log2(1.5))|(get(pval_type) >= pval_thresh)),"black",color)) %>%
                drop_na()
        
                results[[analysis]]<-result
            }
            combined_results <- purrr::reduce(results,function(result1,result2){rbind(result1,result2)})
            print(str(combined_results))    
    
            key_colors <- combined_results$Legend_color
            names(key_colors) <- combined_results$Legend
    
            enhancedplt<-EnhancedVolcano(combined_results,
                            lab =combined_results$mirnaID,
                            x = 'log2FoldChange',
                            y = pval_type, #'padj',#'pvalue',#
                            pCutoff = pval_thresh,
                            FCcutoff = log2(1.5),
                            ylim = c(0, 
                                     ifelse(-log10(pval_thresh)<(max(-log10(combined_results[[pval_type]]))+sd(-log10(combined_results[[pval_type]]))),
                                            max(-log10(combined_results[[pval_type]]))+sd(-log10(combined_results[[pval_type]])),
                                            -log10(pval_thresh))),
                                xlim = c(min(combined_results$log2FoldChange)-sd(combined_results$log2FoldChange),
                                         max(combined_results$log2FoldChange)+sd(combined_results$log2FoldChange)),
                                pointSize = replace((-log10(combined_results[[pval_type]])/max(-log10(combined_results[[pval_type]][!is.infinite(-log10(combined_results[[pval_type]]))])))*5,is.infinite(-log10(combined_results[[pval_type]])),5), #computes ratio between value and max non infinite values times max PointSize, p-values can be 0, thus INf values replaced by max point size,
                            colCustom=key_colors,
                            title =sprintf('Volcano plot with a %s threshold of %s', pval_type, pval_thresh),
                            #subtitle ='NS = non significant',
                           )
            ggsave(filename = file.path(settings[["output_directory"]], "volcano_merged", srna_form, paste0("enhancedCombined_",pval_type,pval_thresh,".pdf")),
                   plot = enhancedplt, height = unit(8, "in"), width = unit(8, "in"))
            }
        }
    }

message("Write Significant result tables according to settings criteria")

write_de_tables <- function(de_results, de_thresholds_settings, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  de_result_thresholds <- purrr::imap(de_results, function(de_result, analysis_name) {
    purrr::imap(de_result, function(result_table, result_name) {
      purrr::imap(de_thresholds_settings, function(criterion, criterion_name) {
        crit <- as.list(criterion)
        x1 <- result_table %>%
          dplyr::filter(padj < !!crit[["padj"]], log2FoldChange > abs(!!crit[["up_log2fc"]])) %>%
          dplyr::mutate(regulation = "upregulated")
        x2 <- result_table %>%
          dplyr::filter(padj < !!crit[["padj"]], log2FoldChange < -abs(!!crit[["down_log2fc"]])) %>%
          dplyr::mutate(regulation = "downregulated")
        x3 <- bind_rows(x1, x2)
        write.csv(
          x = x3,
          file = file.path(
            outdir,
            sprintf("%s-%s-%s.csv", analysis_name, result_name, criterion_name)),
          row.names = FALSE
        )
        x3
      })
    })
  })
  saveRDS(de_result_thresholds, file.path(outdir, "de_result_thresholds.rds"))
  de_result_thresholds
}

de_result_thresholds <- list(
  mature = write_de_tables(
    de_results = de_results$mature,
    de_thresholds_settings = settings[["de_thresholds"]],
    outdir = file.path(settings[["output_directory"]], "de_tables", "mature")
  ),
  hairpin = write_de_tables(
    de_results = de_results$hairpin,
    de_thresholds_settings = settings[["de_thresholds"]],
    outdir = file.path(settings[["output_directory"]], "de_tables", "hairpin")
  )
)


message("Compare gene overlaps with UpSetR")

compare_feat_overlaps <- function(de_result_thresholds, feat_id, gene_overlap_settings, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive=TRUE)
  gene_overlaps <- purrr::map(gene_overlap_settings, function(gene_overlap) {
    name <- gene_overlap[["name"]]
    criterion <- gene_overlap[["criterion"]]
    sets <- gene_overlap[["sets"]]
    gene_sets <- list()
    for (set in sets) {
      set_name <- set[["name"]]
      analysis_name <- set[["set"]][[1]]
      result_name <- set[["set"]][[2]]
      gene_sets[[set_name]] <- de_result_thresholds[[analysis_name]][[result_name]][[criterion]][[feat_id]]
    }

    gene_sets <- purrr::keep(gene_sets, ~length(.) > 0)
    if (length(gene_sets) %in% c(0,1)) {
      saveRDS(gene_sets, file.path(outdir, sprintf("%s-%s_sets.rds", name, feat_id)))
      return(gene_sets)
    }
    pdf(file=file.path(outdir, sprintf("%s.pdf", name)),
        width = 7, height = 4)
    print(
      UpSetR::upset(
        data = UpSetR::fromList(gene_sets),
        nsets = 8, keep.order=TRUE, order.by = "freq"),
      newpage = FALSE
    )
    dev.off()
    saveRDS(gene_sets, file.path(outdir, sprintf("%s-%s_sets.rds", name, feat_id)))
    gene_sets
  })
}

gene_sets = list(
  mature = compare_feat_overlaps(
    de_result_thresholds = de_result_thresholds$mature,
    feat_id = "mirnaID",
    gene_overlap_settings = settings[["gene_overlaps"]],
    outdir =  file.path(settings[["output_directory"]], "gene_overlaps", "mature")
  ),
  hairpin = compare_feat_overlaps(
    de_result_thresholds = de_result_thresholds$hairpin,
    feat_id = "mirnaID",
    gene_overlap_settings = settings[["gene_overlaps"]],
    outdir =  file.path(settings[["output_directory"]], "gene_overlaps", "hairpin")
  )
)



message("Reading miRNA to target genes...")

# Map mirna to target genes:
mirna_target_score <- readr::read_tsv(
  settings["miRDB_prediction"],
  col_names = c("mirnaID", "target", "score"),
  col_types = "ccd"
)

mirna_target_score <- mirna_target_score %>%
  dplyr::filter(score > settings[["mirDB_target_score"]]) %>%
  dplyr::group_by(mirnaID) %>%
  dplyr::mutate(num_targets = n_distinct(target)) %>%
  dplyr::ungroup()

refseq_to_ensembl <- readr::read_csv(
  settings[["refseq_mrna_to_ensembl_gene_id"]],
  col_types = cols()
)

mirna_to_target_mapping <- inner_join(
  mirna_target_score,
  refseq_to_ensembl,
  by = c("target" = "refseq_mrna")
)

# number of targets per mirna:

message("How many mirna have between (0, 10] targets?")
mirna_to_target_mapping %>%
  dplyr::select(mirnaID, num_targets) %>%
  unique() %>%
  dplyr::pull("num_targets") %>%
  cut( breaks = c(0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)) %>%
  table() %>%
  tibble::as_tibble()

mirna_to_target_mapping_filtered <- mirna_to_target_mapping %>%
  dplyr::filter(num_targets < settings[["mirDB_max_targets"]])

if (settings[["organism"]] == "human") {
  gene_symbol_symbol <- rlang::sym("hgnc_symbol")
} else if (settings[["organism"]] == "mouse") {
  gene_symbol_symbol <- rlang::sym("mgi_symbol")
} else {
  stop("Not implemented organism")
}

mirna_to_target_gene_list <- mirna_to_target_mapping_filtered %>%
  group_by(mirnaID) %>%
  summarize(
    ensembl_gene_id = list(ensembl_gene_id),
    gene_symbol = list(!!gene_symbol_symbol),
    num_targets = n(),
    .groups = "drop"
  ) %>%
  ungroup()

de_results_thresholds_with_targets <- purrr::map(
  de_result_thresholds,
  function(de_result_mat_or_hp) {
    purrr::map(
      de_result_mat_or_hp,
      function(de_result_anal) {
        purrr::map(
          de_result_anal,
          function(de_result_crit) {
            purrr::map(
              de_result_crit,
              function(de_result) {
                dplyr::left_join(
                  de_result,
                  mirna_to_target_gene_list,
                  by = "mirnaID"
                )
              }
            )
          }
        )
      }
    )
  }
)
saveRDS(
  de_results_thresholds_with_targets,
  file = file.path(settings[["output_directory"]], "de_results_with_targets.rds")
)


