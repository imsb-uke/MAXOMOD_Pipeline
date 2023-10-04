#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/remove_ensembl_version.R"))
source(here::here("src/shared/R-helpers/create_volcano_plot.R"))
source(here::here("src/shared/R-helpers/memoise_decorator.R"))

suppressPackageStartupMessages({
  library("AnnotationDbi")
  library("biomaRt")
  library("clusterProfiler")
  # library("argparse")
  library("cowplot")
  library("dplyr")
  library("tidyr")
  # library("DESeq2")
  # library("forcats")
  library("ggplot2")
  # library("matrixStats")
  # library("pheatmap")
  library("purrr")
  # library("RColorBrewer")
  library("readr")
  library("tibble")
  library("UpSetR")
  library("EnhancedVolcano")
  library("yaml")
})

settings <- get_settings(omic = "rnaseq", stage = "deg_digestion", dataset_name = NULL) # dataset_name="human")

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


cand_genes <- yaml::read_yaml("candidate-genes.yaml")
cand_genes_species <- cand_genes$sources[[1]]$species
cand_genes_annot_obj <- AnnotationDbi::get(organism_to_annotdbistr(cand_genes_species))
cand_genes_ensembl_src <- clusterProfiler::bitr(
  geneID = cand_genes$sources[[1]]$genes,
  fromType = cand_genes$sources[[1]]$gene_keytype,
  toType = "ENSEMBL",
  OrgDb = cand_genes_annot_obj,
  drop = FALSE
) %>% dplyr::pull("ENSEMBL")


if (settings[["organism"]] == cand_genes_species) {
  cand_genes_ensembl_target <- cand_genes_ensembl_src
} else {
  biomart_convert_ensembl <- function(
    src_species,
    target_species,
    ensembl_src
  ) {
    Sys.sleep(1L)
    src_mart = biomaRt::useMart(
      "ensembl",
      dataset = sprintf("%s_gene_ensembl", cand_genes_species),
      host = "nov2020.archive.ensembl.org"
    )
    target_mart = biomaRt::useMart(
      "ensembl",
      dataset = sprintf("%s_gene_ensembl", settings[["organism"]]),
      host = "nov2020.archive.ensembl.org"
    )
    biomaRt::getLDS(
      attributes = "ensembl_gene_id",
      filters = "ensembl_gene_id",
      values = ensembl_src,
      mart = src_mart,
      attributesL = "ensembl_gene_id",
      martL = target_mart
    ) %>% dplyr::pull(2)
  }
  biomart_convert_ensembl_cache <- memoise_decorator(
    cache_dirname = settings[["cache"]],
    cache_prefix = "biomart_convert_ensembl",
    fun = biomart_convert_ensembl
  )
  cand_genes_ensembl_target <- biomart_convert_ensembl_cache(
    src_species = cand_genes_species,
    target_species = settings[["organism"]],
    ensembl_src = cand_genes_ensembl_src
  )
}



ensg_symbol <- readRDS(file.path(settings[["prefiltpca_dir"]], "ENSG_to_symbol.rds"))

# Helper functions for later using ensg_symbol
# character vector
names_ensg_to_symbol <- function(x) {
  ens <- names(x)
  new_names <- tibble::deframe(ensg_symbol)[ens]
  names(x) <- new_names
  x
}

# data frame with gene_id, will get gene_name column
df_geneid_to_gene_name <- function(x, add = FALSE) {
  if ("gene_name" %in% colnames(x)) {
    return(x)
  }
  out <- left_join(x, ensg_symbol, by = "gene_id")
  if (add) {
    out <- dplyr::select(out, gene_id, gene_name, everything())
  } else {
    out <- dplyr::select(out, gene_name, everything())
  }
  out
}


message("Reading DEG outputs")

de_results <- list()
for (deseq2_setting in settings[["deseq2"]]) {
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

message("Creating single volcano plots")
volcano_dir <- file.path(settings[["output_directory"]], "volcano_single")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
volcano_plots <- purrr::imap(de_results, function(results, deseq2_setting_name) {
  volcano_results_dir <- file.path(volcano_dir, deseq2_setting_name)
  dir.create(volcano_results_dir, showWarnings = FALSE, recursive = TRUE)
  purrr::imap(results, function(result, result_name) {
    result$gene_id_noversion <- remove_ensembl_version(result$gene_id)
    plt <- create_volcano_plot(
      df = result,
      padj_thres = 0.05,
      up_log2fc_thres = log2(1.5),
      down_log2fc_thres = -log2(1.5),
      padj_thres_loose = 0.1,
      label_featurename = "gene_name",
      label_featurefiltername = "gene_id_noversion",
      label_featurefiltervalues = cand_genes_ensembl_target
    )
    ggplot2::ggsave(filename = file.path(volcano_results_dir, sprintf("%s.pdf", result_name)),
           plot = plt, height = unit(15, "in"), width = unit(15, "in"))
    ggplot2::ggsave(filename = file.path(volcano_results_dir, sprintf("%s-small.pdf", result_name)),
           plot = plt, height = unit(5, "in"), width = unit(7, "in"))
    plt
  })
})

message("Creating merged volcano plots")
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

message("Creating merged enhanced volcano plots")
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
                        lab =combined_results$gene_name,
                        x = 'log2FoldChange',
                        y = pval_type,
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
                        title = sprintf('Volcano plot with a %s threshold of %s', pval_type , pval_thresh),
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
      x3 <- dplyr::bind_rows(x1, x2)
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


message("Compare gene overlaps with UpSetR")

gene_overlaps_dir <- file.path(settings[["output_directory"]], "gene_overlaps")
dir.create(gene_overlaps_dir, showWarnings = FALSE, recursive=TRUE)
gene_overlaps <- purrr::map(settings[["gene_overlaps"]], function(gene_overlap) {
  name <- gene_overlap[["name"]]
  criterion <- gene_overlap[["criterion"]]
  sets <- gene_overlap[["sets"]]
  gene_sets <- list()
  for (set in sets) {
    set_name <- set[["name"]]
    analysis_name <- set[["set"]][[1]]
    result_name <- set[["set"]][[2]]
    tmp <- de_result_thresholds[[analysis_name]][[result_name]][[criterion]][["gene_id"]]
    if (length(tmp) > 0){
        gene_sets[[set_name]] <- tmp
        }
  }
  if (length(gene_sets) > 1) {
  pdf(file=file.path(gene_overlaps_dir, sprintf("%s.pdf", name)),
      width = 7, height = 4)
  print(
    UpSetR::upset(
      data = UpSetR::fromList(gene_sets),
      nsets = 8, keep.order=TRUE, order.by = "freq"),
    newpage = FALSE
  )
  dev.off()
    }
  saveRDS(gene_sets, file.path(gene_overlaps_dir, sprintf("%s-gene_sets.rds", name)))
  gene_sets
})
