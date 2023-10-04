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

make_go_enr = function(genelist, prefix, folder, orgDB, pvalueCutoff, qvalueCutoff,
    pAdjustMethod) {
    ego <- enrichGO(gene = genelist$ENTREZID, OrgDb = orgDB, ont = prefix, pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = T)

    if (!is.null(ego) && nrow(ego) != 0) {
        p1 = dotplot(ego, showCategory = 30)
        ggsave(p1, filename = file.path(folder, paste0("dotplot_", prefix, ".png")),
            width = 20, units = "cm")

        p2 = barplot(ego, showCategory = 30)
        ggsave(p2, filename = file.path(folder, paste0("barplot_", prefix, ".png")),
            width = 20, units = "cm")
        write.csv(data.frame(ego), file = file.path(folder, paste0("pathways_", prefix,
            ".csv")))
    } else {
        file.create(file.path(folder, paste0("pathways_", prefix, ".csv")))
    }
}

make_kegg_enr = function(genelist, folder, organism, pvalueCutoff, qvalueCutoff,
    pAdjustMethod) {
    ego <- enrichKEGG(gene = genelist$ENTREZID, organism = organism, qvalueCutoff = qvalueCutoff,
        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod)

    if (!is.null(ego) && nrow(ego) != 0) {
        p1 = dotplot(ego, showCategory = 30)
        ggsave(p1, filename = file.path(folder, "dotplot.png"), width = 20, units = "cm")

        p2 = barplot(ego, showCategory = 30)
        ggsave(p2, filename = file.path(folder, "barplot.png"), width = 20, units = "cm")
        write.csv(data.frame(ego), file = file.path(folder, "pathways.csv"))
    } else {
        file.create(file.path(folder, "pathways.csv"))
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

print(paste0(dim(result_table)[1], " genes left after mapping ENTREZ Ids"))

if (!("SYMBOL" %in% colnames(result_table))) {
    mapping <- mapIds(orgDB_set, keys = result_table$ENTREZID, column = c("SYMBOL"),
        keytype = "ENTREZID", multiVals = "first")
    mapping = unlist(mapping, use.names = TRUE)
    mapping = data.frame(ENTREZID = names(mapping), SYMBOL = mapping)
    result_table <- merge(result_table, mapping, by = c("ENTREZID"))
}

result_table = result_table[!duplicated(result_table[, "SYMBOL"]), ]

print(paste0(dim(result_table)[1], " genes left after removing duplicates"))


rownames(result_table) = result_table$SYMBOL


# Filter

genelist = filter(result_table, eval(parse(text = params[["filter"]])))

print(paste0(dim(genelist)[1], " genes left after filtering"))

# Make GO enrichment plots

gopath = file.path(output, "go")
dir.create(gopath, recursive = TRUE, showWarnings = FALSE)

make_go_enr(genelist, "BP", gopath, orgDB_set, pvalueCutoff = params[["go_pvalueCutoff"]],
    qvalueCutoff = params[["go_qvalueCutoff"]], pAdjustMethod = params[["go_pAdjustMethod"]])
make_go_enr(genelist, "MF", gopath, orgDB_set, pvalueCutoff = params[["go_pvalueCutoff"]],
    qvalueCutoff = params[["go_qvalueCutoff"]], pAdjustMethod = params[["go_pAdjustMethod"]])
make_go_enr(genelist, "CC", gopath, orgDB_set, pvalueCutoff = params[["go_pvalueCutoff"]],
    qvalueCutoff = params[["go_qvalueCutoff"]], pAdjustMethod = params[["go_pAdjustMethod"]])

print("GO enrichment done")

# Make KEGG enrichment plots

keggpath = file.path(output, "kegg")
dir.create(keggpath, recursive = TRUE, showWarnings = FALSE)

ego = make_kegg_enr(genelist, keggpath, organism = params[["stages"]][[dataset_name]][["type"]],
    pvalueCutoff = params[["kegg_pvalueCutoff"]], qvalueCutoff = params[["kegg_qvalueCutoff"]],
    pAdjustMethod = params[["kegg_pAdjustMethod"]])

setwd(keggpath)

for (i in 1:length(ego$ID)) {
    pv.out <- pathview(gene.data = result_table[params[["color"]]], pathway.id = ego$ID[i],
        species = params[["stages"]][[dataset_name]][["type"]], out.suffix = dataset_name,
        kegg.native = T, gene.idtype = "SYMBOL", low = "blue", mid = "gray", high = "red",
        same.layer = F)

    file.remove(paste0(ego$ID[i], ".png"))
    file.remove(paste0(ego$ID[i], ".xml"))
}

print("kegg enrichment done")
