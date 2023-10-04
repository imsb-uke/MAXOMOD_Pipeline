library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(tidyverse)
library(argparse)
source(file.path("src","shared","WGCNA","perform_WGCNA_functions.R"))

get_args <- function() {

    # create parser object
    parser <- ArgumentParser()
    parser$add_argument("--dataset-name", type = "character", default = NULL, help = "dataset key from params.yaml to use")


    args <- parser$parse_args()
    args
}

# MAIN


## Read params

args = get_args()

dataset_name = args$dataset_name[1]
print(dataset_name)

params = yaml::read_yaml("params.yaml")[["WGCNA"]][["proteomics"]][["stages"]][[dataset_name]]



if (params[["organism"]]=="Homo sapiens"){
    OrgDb<-org.Hs.eg.db
} else if(params[["organism"]]=="Mus musculus"){
    OrgDb<-org.Mm.eg.db
} else {stop(paste0("Not a valid organism! ",params[["organism"]]))}

vst_mat = read.csv(params[["input_file"]])
rownames(vst_mat) = vst_mat$X
vst_mat<-vst_mat[2:length(vst_mat)]

#print(dim(vst_mat))
#print(colnames(vst_mat))
#print(head(vst_mat))

#Read in protein to ensembl id mapping
prot_to_gene <- readr::read_tsv(params[["uniprot_to_ensembl"]], col_names = FALSE)
prot_to_gene <- prot_to_gene[,c("X2", "X19")]
colnames(prot_to_gene) <- c("UniProtName", "gene_id")
prot_to_gene<-prot_to_gene[!is.na(prot_to_gene)[,c("gene_id")],]
#print(dim(prot_to_gene))

#get gene symbol from ensembl id
mapping <- mapIds(OrgDb, keys = pull(prot_to_gene,"gene_id"), column = c("SYMBOL"), keytype = c("ENSEMBL"), multiVals="first")
prot_to_gene<-cbind(prot_to_gene,data.frame("gene_name"=mapping))

prot_to_gene = prot_to_gene[!duplicated(prot_to_gene$gene_name),]

#intersect selects all symbols of genes that were already prefiltered in expression data
matching_pairs = intersect(colnames(vst_mat), prot_to_gene$UniProtName)
#print(matching_pairs)
#print(length(matching_pairs))
rownames(prot_to_gene) = prot_to_gene$UniProtName
prot_to_gene_sub = prot_to_gene[matching_pairs,]
vst_mat = t(vst_mat)[matching_pairs,]

#print(dim(prot_to_gene_sub))
#print(dim(vst_mat))
#print(head(vst_mat))

#Checking if all rownames are matching, then replace ensemble gene id with gene_name symbols
if (all(rownames(vst_mat) == prot_to_gene_sub$UniProtName)){
  rownames(vst_mat) = prot_to_gene_sub$gene_name 
} else {
  stop("Rownames are not matching!")
}

anno = read.csv(params[["annotations"]])
#print(head(anno))
#print(head(vst_mat))
#print(dim(vst_mat))

stopifnot(all((colnames(vst_mat)) == (anno$SampleID)))

#creating conditions
conditions<-factor(sapply(rownames(anno),function(id){paste(anno[id,c("Condition","Sex")],collapse = '_')}))
#create trait dataframe, where each sample is a vector of conditions where the trait is set to 1, the other traits to 0
trait_df<-data.frame(t(sapply(conditions,function(condition){table(condition)})),row.names=anno$SampleID)

#print(trait_df)

#create results folder
out_dir<-params[["output_dir"]]
de_dir<-params[["de_dir"]]

#load and format text mining literature results
literature_path<-yaml::read_yaml("params.yaml")[["WGCNA"]][["proteomics"]][["literature_dir"]]
literatureDF<-read.csv(literature_path,sep='|')
literature_entrez_ids<-sapply(literatureDF[,"term"],function(term){strsplit(term, split = ":")[[1]][2]})
literatureDF["term"]<-mapIds(OrgDb, keys=literature_entrez_ids, column=c("SYMBOL"), 
                            keytype="ENTREZID", multiVals="first")

#load als markers
als_markers_source <- yaml::read_yaml("params.yaml")[["WGCNA"]][["proteomics"]][["als_markers_source"]]
als_markers <- read.csv(als_markers_source) %>%
    distinct(marker_genes) %>%
    select(marker_genes)

perform_WGCNA(vst_mat,trait_df,out_dir,de_dir,literatureDF,als_markers,mapping=prot_to_gene_sub,OrgDb=OrgDb)