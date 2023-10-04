#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/create_volcano_plot.R"))

suppressPackageStartupMessages({
  library("cowplot")
  library("dplyr")
  library("ggplot2")
  library("purrr")
  library("UpSetR")
  library("EnhancedVolcano")
  library("stringr")
})

settings <- get_settings(omic = "phosphoproteomics", stage = "deg_digestion", dataset_name = NULL)

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

message("Loading differential expression of proteins...")
de_results <- readRDS(settings[["de_results_rds"]])

message("Creating volcano plots")
volcano_dir <- file.path(settings[["output_directory"]], "volcano_single")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
volcano_plots <- purrr::imap(de_results, function(results, analysis_name) {
  volcano_results_dir <- file.path(volcano_dir, analysis_name)
  dir.create(volcano_results_dir, showWarnings = FALSE, recursive = TRUE)
  purrr::imap(results, function(result, result_name) {
    plt <- create_volcano_plot(
      df = result,
      padj_thres = 0.05,
      up_log2fc_thres = log2(1.2),
      down_log2fc_thres = -log2(1.2),
      padj_thres_loose = 0.1
    )
    ggsave(filename = file.path(volcano_results_dir, sprintf("%s.pdf", result_name)),
           plot = plt, height = unit(4, "in"), width = unit(6, "in"))
    #create additional enhancedVolcano plots
    for (pval_thresh in c(0.05,0.1)) {
    enhancedplt<-EnhancedVolcano(result,
                    lab =result$UniProtName,
                    x = 'log2FoldChange',
                    y = 'padj',#'pvalue',#
                    pCutoff = pval_thresh,
                    FCcutoff = log2(1.2),
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
    #return first images
    plt
  })
})

volcano_dir <- file.path(settings[["output_directory"]], "volcano_merged")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
for (volcano_merged_plot in settings[["volcano_plots"]]) {
  name <- volcano_merged_plot[["name"]]
  subplots <- volcano_merged_plot[["subplots"]]
  plot_list <- list()
  for (subplot in subplots) {
    analysis <- subplot[["name"]][[1]]
    result <- subplot[["name"]][[2]]
    subplot_title <- subplot[["title"]]
    plt <- volcano_plots[[analysis]][[result]]
    plt <- plt + labs(title = subplot_title)
    plot_list <- append(plot_list, list(plt))
  }
  plt <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  ggsave(filename = file.path(volcano_dir, sprintf("%s.pdf", name)),
         plot = plt, height = unit(6, "in"), width = unit(8, "in"))
  plt
}

#Merge male and female volcano plot into one
for (pval_type in c('padj','pvalue')){
    for (pval_thresh in c(0.05,0.1)) { 
        results <- list()
        for (setting in settings[["merged_volcano"]]){
            analysis <- setting[["analysis"]]
            title <- setting[["title"]]
            result_name <- setting[["result_name"]]
            color <- setting[["color"]]
            result <- de_results[[analysis]][[result_name]] %>%
                mutate(Legend=if_else((between(log2FoldChange, -log2(1.2), log2(1.2))|(get(pval_type) >= pval_thresh)),"not significant",title),
                      Legend_color=if_else((between(log2FoldChange, -log2(1.2), log2(1.2))|(get(pval_type) >= pval_thresh)),"black",color))
        
            results[[analysis]]<-result
        }
        combined_results <- purrr::reduce(results,function(result1,result2){rbind(result1,result2)})
        print(str(combined_results))    
    
        key_colors <- combined_results$Legend_color
        key_colors <- combined_results$Legend_color
        names(key_colors) <- combined_results$Legend
    
        enhancedplt<-EnhancedVolcano(combined_results,
                        lab =combined_results$UniProtName,
                        x = 'log2FoldChange',
                        y = pval_type,
                        pCutoff = pval_thresh,
                        FCcutoff = log2(1.2),
                        ylim = c(0, 
                                 ifelse(-log10(pval_thresh)<(max(-log10(combined_results[[pval_type]]))+sd(-log10(combined_results[[pval_type]]))),
                                    max(-log10(combined_results[[pval_type]]))+sd(-log10(combined_results[[pval_type]])),
                                    -log10(pval_thresh))),
                        xlim = c(min(combined_results$log2FoldChange)-sd(combined_results$log2FoldChange),
                                 max(combined_results$log2FoldChange)+sd(combined_results$log2FoldChange)),
                        pointSize = replace((-log10(combined_results[[pval_type]])/max(-log10(combined_results[[pval_type]][!is.infinite(-log10(combined_results[[pval_type]]))])))*5,is.infinite(-log10(combined_results[[pval_type]])),5), #computes ratio between value and max non infinite values times max PointSize, p-values can be 0, thus INf values replaced by max point size,
                        colCustom=key_colors,
                        title = sprintf('Volcano plot with a %s threshold of %s', pval_type , pval_thresh),,
                        #subtitle ='NS = non significant',
                       )
        ggsave(filename = file.path(volcano_dir, paste0("enhancedCombined_",pval_type,pval_thresh,".pdf")),
               plot = enhancedplt, height = unit(8, "in"), width = unit(8, "in"))
        }
    }

message("Write Significant result tables according to settings criteria")

de_tables_outdir <- file.path(settings[["output_directory"]], "de_tables")
dir.create(de_tables_outdir, showWarnings = FALSE, recursive = TRUE)
de_result_thresholds <- purrr::imap(de_results, function(de_result, analysis_name) {
  purrr::imap(de_result, function(result_table, result_name) {
    purrr::imap(settings[["de_thresholds"]], function(criterion, criterion_name) {
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
          de_tables_outdir,
          sprintf("%s-%s-%s.csv", analysis_name, result_name, criterion_name)),
        row.names = FALSE
      )
      x3
    })
  })
})
saveRDS(de_result_thresholds, file.path(de_tables_outdir, "de_result_thresholds.rds"))


message("Generating UpSetR plots of protein sets...")


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
  gene_overlaps
}

gene_sets <- compare_feat_overlaps(
  de_result_thresholds = de_result_thresholds,
  feat_id = "UniProtName",
  gene_overlap_settings = settings[["gene_overlaps"]],
  outdir =  file.path(settings[["output_directory"]], "gene_overlaps")
)

message("DONE")
