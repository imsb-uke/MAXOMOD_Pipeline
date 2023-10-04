suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(DOSE)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(pathview)
    library(tidyverse)
})

make_kegg_gsea = function(genelist, folder, organism, pvalueCutoff, pAdjustMethod,
    scoreType) {
    ego <- gseKEGG(geneList = genelist, organism = organism, pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod, scoreType = scoreType, keyType = "ncbi-geneid")

    if (!is.null(ego) && nrow(ego) != 0) {
        p1 = dotplot(ego, showCategory = 30)
        ggsave(p1, filename = file.path(folder, "dotplot.png"), width = 20, units = "cm")
        write.csv(data.frame(ego), file = file.path(folder, "pathways.csv"))
    }
    ego
}


get_args <- function() {

    # create parser object
    parser <- argparse::ArgumentParser()
    parser$add_argument("--category", type = "character", default = NULL, help = "Category key from params.yaml to use")
    parser$add_argument("--dataset-name", type = "character", default = NULL, help = "dataset key from params.yaml to use")
    parser$add_argument("--omic", type = "character", default = NULL, help = "omic key from params.yaml to use")
    parser$add_argument("--method", type = "character", default = NULL, help = "method key from params.yaml to use")


    args <- parser$parse_args()
    args
}

# MAIN


## Read params

args = get_args()

category = args$category[1]
method = args$method[1]
omic = args$omic[1]
dataset_name = args$dataset_name[1]


params = yaml::read_yaml("params.yaml")[[category]][[omic]][[method]]


orgDB_sets = c(mmu = org.Mm.eg.db, hsa = org.Hs.eg.db)
taxid = c(mmu = 10090, hsa = 9606)



infile = params[["stages"]][[dataset_name]][["infile"]]
output = params[["stages"]][[dataset_name]][["output"]]

orgDB_set = orgDB_sets[[params[["stages"]][[dataset_name]][["type"]]]]

# Read data

result_table = read.csv(infile)
print(paste0("Loaded ", dim(result_table)[1], " genes"))


# Map IDs

genecol = params[["genecol"]]
keytype = params[["stages"]][[dataset_name]][["keytype"]]

if (keytype == "ENSEMBL") {
    result_table[, genecol] <- gsub("\\.[0-9]*$", "", result_table[, genecol])
}

if (params[["idtype"]] == "protein") {
    mapping = read.csv(params[["stages"]][[dataset_name]][["protein_id_mapping"]],
        stringsAsFactors = FALSE, sep = "\t", header = FALSE)[, c(1, 2)]
    colnames(mapping) = c("UNIPROT", genecol)
    result_table <- merge(result_table, mapping, by = c(genecol))
    genecol = "UNIPROT"
    keytype = "UNIPROT"
    print(paste0(dim(result_table)[1], " genes left after protein id mapping"))
}

mapping <- mapIds(orgDB_set, keys = result_table[, genecol], column = c("ENTREZID"),
    keytype = keytype, multiVals = params[["multiVals"]])
mapping = mapping[lapply(mapping, length) > 0]
mapping = data.frame(mapping)
colnames(mapping) = c("ENTREZID")
mapping[, genecol] = rownames(mapping)
result_table <- merge(result_table, mapping, by = c(genecol))

print(paste0(dim(result_table)[1], " genes left after mapping ENTREZIDs"))

result_table = result_table[!duplicated(result_table[, "ENTREZID"]), ]
result_table = result_table[!is.na(result_table[, "ENTREZID"]), ]

print(paste0(dim(result_table)[1], " genes left after removing duplicates"))


# Filter

genelist = result_table[[params[["column"]]]]
names(genelist) = result_table$ENTREZID
genelist = sort(genelist, decreasing = T)


# Make GSEA enrichment plots

gseapath = file.path(output, "gsea")
dir.create(gseapath, recursive = TRUE, showWarnings = FALSE)

ego = make_kegg_gsea(genelist, gseapath, organism = params[["stages"]][[dataset_name]][["type"]],
    pvalueCutoff = params[["pvalueCutoff"]], pAdjustMethod = params[["pAdjustMethod"]],
    scoreType = params[["scoreType"]])

print("GSEA enrichment done")
