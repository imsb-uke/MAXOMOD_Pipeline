#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

library("GeneTonic")
library("topGO")
library(DESeq2)
library(pcaExplorer)
library(dplyr)
library("ggplot2")

run_gene_tonic = function(de_tab, female_dds, cat, fname, organism, foldchange, padj) {
    rownames(de_tab) = de_tab$gene_id
    
    female_dds = female_dds[de_tab$gene_id, ]
    
    all((rownames(female_dds)) == (de_tab$gene_id))
    
    de_tab_filt = de_tab[abs(de_tab$log2FoldChange) > foldchange & de_tab$padj < 
        padj, ]
    de_tab_filt = as.data.frame(de_tab_filt[!is.na(de_tab_filt$padj), ])
    
    
    bg_ids <- de_tab$gene_name[rowSums(counts(female_dds)) > 0]
    
    topgoDE_female <- pcaExplorer::topGOtable(de_tab_filt$gene_name, bg_ids, ontology = cat, 
        mapping = organism, geneID = "symbol", topTablerows = 500)
    
    res_enrich_female <- shake_topGOtableResult(topgoDE_female)
    
    anno_df <- data.frame(gene_id = de_tab$gene_id, gene_name = de_tab$gene_name, 
        stringsAsFactors = FALSE, row.names = de_tab$gene_id)
    
    res_enrich_female <- get_aggrscores(res_enrich = res_enrich_female, res_de = de_tab_filt, 
        annotation_obj = anno_df, aggrfun = mean)
    
    
    pdf(fname)
    print(gs_volcano(res_enrich_female, p_threshold = padj, color_by = "aggr_score", 
        volcano_labels = 10, gs_ids = NULL, plot_title = paste0("RNA Seq ", cat)))
    dev.off()
    
    
    return(res_enrich_female)
}

gs_summary_overview_pair_zscore = function(res_enrich, res_enrich2, n_gs = 50, p_value_column = "gs_pvalue", 
    color_by = "z_score", alpha_set2 = 1) {
    gs_set1 <- res_enrich$gs_id
    gs_set2 <- res_enrich2$gs_id
    gs_common <- intersect(gs_set1, gs_set2)
    if (length(gs_common) == 0) {
        stop("No gene sets have been found in common to the two enrichment results")
    }
    
    gs_common <- gs_common[seq_len(min(n_gs, length(gs_common)))]
    common_re1 <- res_enrich[gs_common, ]
    common_re2 <- res_enrich2[gs_common, ]
    
    common_re1[[p_value_column]][is.na(common_re1[[p_value_column]])] = 1
    common_re2[[p_value_column]][is.na(common_re2[[p_value_column]])] = 1
    
    common_re1$minus_logp10 <- -log10(common_re1[[p_value_column]])
    common_re2$minus_logp10 <- -log10(common_re2[[p_value_column]])
    re_both <- common_re1
    re_both[["minus_logp10_2"]] <- common_re2$minus_logp10
    re_both[[color_by]] <- common_re1[[color_by]]
    re_both[[paste0(color_by, "_2")]] <- common_re2[[color_by]]
    re_both[["zscore_diff"]] = abs(common_re1[[color_by]] - common_re2[[color_by]])
    
    re_both_sorted <- re_both %>% arrange(.data$zscore_diff) %>% mutate(gs_description = factor(.data$gs_description, 
        .data$gs_description))
    
    p <- ggplot(re_both_sorted, aes_string(x = "gs_description", y = "z_score")) + 
        geom_segment(aes_string(x = "gs_description", xend = "gs_description", y = "z_score_2", 
            yend = "z_score"), color = "grey") + geom_point(aes(fill = .data[["minus_logp10"]]), 
        size = 4, pch = 21) + geom_point(aes_string(y = "z_score_2", col = "minus_logp10_2"), 
        size = 4, alpha = alpha_set2) + scale_color_gradient2(low = "#313695", mid = "#FFFFE5", 
        high = "#A50026") + scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", 
        high = "#A50026", guide = FALSE) + coord_flip() + labs(x = "Gene set description", 
        y = "Z-Score", col = "logp10") + ylim(0, NA) + theme_minimal()
    
    return(p)
}


##### MAIN #####

settings <- get_settings(omic = "rnaseq", stage = "genetonic", dataset_name = NULL)  #dataset_name = 'c9')

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

if (settings[["organism"]] == "mmusculus") {
    organism <- "org.Mm.eg.db"
} else if (settings[["organism"]] == "hsapiens") {
    organism <- "org.Hs.eg.db"
} else {
    stop("Not implemented error")
}


categories = settings[["categories"]]  #c('MF','BP','CC')

annotation = read.table(settings[["annotation"]], sep = ",", header = T, row.names = 1)

results = list()

for (sex in c("male", "female")) {
    
    results[[sex]] = list()
    dds = readRDS(file.path(settings[[paste0("inputs_", sex)]], "dds.rds"))
    de_tab = read.csv(file.path(settings[[paste0("inputs_", sex)]], "mut_vs_ctrl.csv"))
    
    for (cat in categories) {
        fname = file.path(settings[["output_directory"]], paste0(sex, "_de_topGO_", 
            cat, "_volcano.pdf"))
        res = run_gene_tonic(de_tab, dds, cat, fname, organism, foldchange = settings[["foldchange_threshold"]], 
            padj = settings[["padj_threshold"]])
        write.csv(res, file = file.path(settings[["output_directory"]], paste0(sex, 
            "_de_topGO_", cat, ".csv")))
        results[[sex]][[cat]] = res
    }
}


for (cat in categories) {
    fname = file.path(settings[["output_directory"]], paste0("comparison_enrichment_", 
        cat, ".pdf"))
    
    pdf(fname, width = 20)
    
    print(gs_summary_overview_pair(res_enrich = results[["female"]][[cat]], res_enrich2 = results[["male"]][[cat]], 
        n_gs = settings[["nterms"]]))
    
    dev.off()
    fname = file.path(settings[["output_directory"]], paste0("comparison_enrichment_", 
        cat, "_zscore.pdf"))
    
    
    pdf(fname, width = 20)
    
    print(gs_summary_overview_pair_zscore(res_enrich = results[["female"]][[cat]], 
        res_enrich2 = results[["male"]][[cat]], n_gs = settings[["nterms"]]))
    
    dev.off()
}
