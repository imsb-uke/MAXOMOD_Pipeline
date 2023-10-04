library(tidyverse)
library(jsonlite)
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

library(argparse)

library(WGCNA)
library(flashClust)

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

params = yaml::read_yaml("params.yaml")[["WGCNA"]][["rnaseq"]][["stages"]][[dataset_name]]



if (params[["organism"]]=="Homo sapiens"){
    OrgDb<-org.Hs.eg.db
    organism = "hsa"
} else if(params[["organism"]]=="Mus musculus"){
    OrgDb<-org.Mm.eg.db
    organism = "mmu"
} else {stop(paste0("Not a valid organism! ",params[["organism"]]))}



# Set number of threads
nThreads = 30
vst_mat = read.csv(params[["input_file"]])
rownames(vst_mat) = vst_mat$X
vst_mat<-vst_mat[2:length(vst_mat)]

dir.create(params[["output_dir"]], showWarnings =FALSE, recursive = TRUE)

ensg_symbol = read.csv(params[["ensembl_to_symbol"]])

#Replace gene identifiers from ensemble id to symbol gene_name
ensg_symbol = ensg_symbol[!duplicated(ensg_symbol$gene_name),]

#intersect selects all symbols of genes that were already prefiltered in expression data
matching_pairs = intersect(rownames(vst_mat), ensg_symbol$gene_id)
rownames(ensg_symbol) = ensg_symbol$gene_id
ensg_symbol_sub = ensg_symbol[matching_pairs,]
vst_mat = vst_mat[matching_pairs,]

#Checking if all rownames are matching, then replace ensemble gene id with gene_name symbols
if (all(rownames(vst_mat) == ensg_symbol_sub$gene_id)){
  rownames(vst_mat) = ensg_symbol_sub$gene_name
} else {
  stop("Rownames are not matching!")
}

#Creating annotation object -> 4 conditions (ctrl_male,ctrl,female,mut_male,mut_female)
anno = read.csv(params[["annotations"]],
               row.names=1)
colnames(vst_mat) <- gsub("\\.", "-", colnames(vst_mat))

stopifnot(all((colnames(vst_mat)) == (anno$SampleID)))

#creating conditions
conditions<-factor(sapply(rownames(anno),function(id){paste(anno[id,c("Condition","Sex")],collapse = '_')}))
#create trait dataframe, where each sample is a vector of conditions where the trait is set to 1, the other traits to 0
trait_df<-data.frame(t(sapply(conditions,function(condition){table(condition)})),row.names=anno$SampleID)


# This creates an object called "datExpr" that contains the normalized counts file output from DESeq2
datExpr = vst_mat
# "head" the file to preview it
head(datExpr) # You see that genes are listed in a column named "X" and samples are in columns

# Manipulate file so it matches the format WGCNA needs
#row.names(datExpr) = datExpr$X
#datExpr$X = NULL
datExpr = t(datExpr) # now samples are rows and genes are columns

dim(datExpr) # 48 samples and 1000 genes (you will have many more genes in reality)
head(datExpr)
# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

#if (!gsg$allOK)
#	{if (sum(!gsg$goodGenes)>0)
#		printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
#		if (sum(!gsg$goodSamples)>0)
#			printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
#		datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
#		}

#Create an object called "datTraits" that contains your trait data


#datTraits = anno[,-1]
#datTraits = datTraits[,c("control","trait")]
#head(datTraits)
datTraits = trait_df

#form a data frame analogous to expression data that will hold the clinical traits.
#rownames(datTraits) = datTraits$nf_id
#datTraits$nf_id = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

# Cluster samples by expression ----------------------------------------------------------------

A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors((datTraits),signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

png(filename=file.path(params[["output_dir"]],"sampleDendorgramAndTraitHeatmap.png"))
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")
dev.off()

# Choose a soft threshold power
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9

png(file.path(params[["output_dir"]],"scaleIndependence.png"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
dev.off()

png(file.path(params[["output_dir"]],"meanConnectivity.png"))
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
# ----------------------------------------------------------------------------
#from this plot, we would choose a power of 10 becuase it's the lowest power for which the scale free topology index reaches 0.90

enableWGCNAThreads(nThreads=nThreads)


softPower = params[["softPower"]]
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
arrow::write_feather(as.data.frame(adjacency), file.path(params[["output_dir"]],"adjacency_matrix.feather"))

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# Generate Modules --------------------------------------------------------

# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

png(file.path(params[["output_dir"]],"geneClusteringTOM.png"))
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()

#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors)#,softPower = 14)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= file.path(params[["output_dir"]],"Network_allSamples_signed_RLDfiltered.RData"))

#plots tree showing how the eigengenes cluster together
png(file.path(params[["output_dir"]],paste("eigenGenesClustering_", 6, ".png", sep="")))
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()

#softPower = 10 # see scale independence file "softthreshold_SOD1.png"



#### -----------------------------------------------------------------------
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = params[["MEDissThres"]]
print(MEDissThres)
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
png(file.path(params[["output_dir"]],"geneTree.png"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= file.path(params[["output_dir"]],"Network_allSamples_signed_nomerge_RLDfiltered.RData"))


# Correlate traits --------------------------------------------------------

#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use="p")

write.csv(moduleTraitCor, file.path(params[["output_dir"]],"moduleTraitCor.csv"))

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

write.csv(moduleTraitPvalue, file.path(params[["output_dir"]],"moduleTraitPvalue.csv"))

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), " (",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar=  c(6, 3, 3, 3))
png(file.path(params[["output_dir"]],"moduleTraitRelationshipHeatmap_all.png"))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               #ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               #textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
dev.off()
#gets modules where the mut conditions are highly correleted

moduleSelectionThreshold = params[["moduleSelectionThreshold"]]

module_selection<-(moduleTraitCor[,params[["selection"]][1]]<moduleSelectionThreshold&moduleTraitCor[,params[["selection"]][2]]<moduleSelectionThreshold)&(moduleTraitCor[,params[["selection"]][3]]>-moduleSelectionThreshold|moduleTraitCor[,params[["selection"]][4]]>-moduleSelectionThreshold)
#display the corelation values with a heatmap plot
par(mar=  c(6, 9, 3, 3))
png(file.path(params[["output_dir"]],"moduleTraitRelationshipHeatmap_corrAll.png"), width=600, height=600)
labeledHeatmap(Matrix= t(moduleTraitCor[module_selection,]),
               yLabels= names(datTraits),
               xLabels= names(MEs)[module_selection],
               xSymbols= names(MEs)[module_selection],
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               #textMatrix= textMatrix[module_selection,],
               setStdMargins= TRUE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
dev.off()

#get selected genes
selected_Colors<-sapply(names(MEs)[module_selection],function(ME){sub("ME","",ME)})
#names(datExpr)[moduleColors=="brown"]
module_genes<-sapply(selected_Colors,function(module_color){rownames(t(datExpr)[mergedColors==module_color,])})

saveRDS(module_genes,file=file.path(params[["output_dir"]],'module_genes.rds'))
write(toJSON(module_genes),file=file.path(params[["output_dir"]],'module_genes.json'))

## enrichment of each module



make_go_enr = function(genelist,
                       fname,
                       folder,
                       OrgDb,
                       fromType="SYMBOL",
                       prefixes=c("BP","MF","CC","All"),
                       image_width=30){

  gene.df <- bitr((genelist), fromType = fromType,
                  toType = c("ENTREZID"),
                  OrgDb = OrgDb)
  for (prefix in prefixes) {
    ego <- enrichGO(gene          = gene.df$ENTREZID,
                    OrgDb         = OrgDb,
                    ont           = prefix,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable = T)

    dir.create(file.path(folder, fname))

    if (!is.null(ego) && nrow(ego)!=0){
      write.csv(ego, file.path(folder, fname, paste0("go_enr_", prefix,"_",fname,".csv")))

      p1 = dotplot(ego, showCategory=10) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p1, filename = file.path(folder,fname, paste0("dotplot_go_", prefix,"_",fname,".jpg")), width = image_width, units="cm")

      p2 = barplot(ego, showCategory=10) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p2, filename = file.path(folder, fname, paste0("barplot_go_", prefix,"_",fname,".jpg")), width = image_width, units="cm")
    }
  }
}

make_kegg_enr <- function(deg_input,
                          fname,
                          folder,
                          OrgDb,
                          organism,
                          pvalueCutoff=0.05,
                          FC_col="log2FoldChange",
                          degs=NULL){

  old_wd <- getwd()
  dir.create(file.path(folder, 'kegg'))

  if (length(row.names(deg_input)) == 0) {
        return(NULL)
        }

  genes_entrezid <-  mapIds(OrgDb, keys=row.names(deg_input), column=c("ENTREZID"),
                            keytype="SYMBOL", multiVals="first")

  ego <- enrichKEGG(gene=genes_entrezid,
                    organism = organism,
                    pvalueCutoff  = pvalueCutoff)

  write.csv(ego, file.path(folder, 'kegg', paste0("kegg_enr_",fname,".csv")))


  if (!is.null(degs)){
    deg_input<-rbind(deg_input,degs)
  }

  if (!is.null(ego) && nrow(ego)!=0){
    setwd(file.path(folder, 'kegg'))
    for(i in 1:length(ego$ID))
    {
      pv.out <- pathview(
        gene.data  = deg_input[FC_col],
        pathway.id = ego$ID[i],
        species = organism,
        out.suffix = fname,#ego$ID[i],
        kegg.native = T,
        gene.idtype="SYMBOL",
        low = "blue",
        mid = "gray",
        high = "red",
        same.layer = F)
    }
    setwd(old_wd)
  }
}


overview_df<- data.frame()


for (module in names(module_genes)){
  make_go_enr(module_genes[[module]],module,params[["output_dir"]],OrgDb=OrgDb)

  ## kegg enrichment
  de_symbol_col<-"gene_name"
  types = c("only_females", "only_males")[c(moduleTraitCor[module,params[["selection"]][3]]>-moduleSelectionThreshold,moduleTraitCor[module,params[["selection"]][4]]>-moduleSelectionThreshold)]
  for(ll in types){
    obj <- file.path(params[["de_dir"]],ll,"mut_vs_ctrl.csv")
    de_results <- read.csv(obj)
    de_results<-na.omit(de_results[ !duplicated(de_results[,de_symbol_col]), ])
    rownames(de_results) = de_results[,de_symbol_col]

    make_kegg_enr(de_results[rownames(de_results) %in% module_genes[[module]],],
                  paste0(module,'_',ll),
                  OrgDb=OrgDb,
                  organism=organism,
                  file.path(params[["output_dir"]],module),
                  degs=na.omit(de_results[de_results["padj"]<=0.05,]))
  }
  ## get module overview table
  num_module_genes <- length(module_genes[[module]])
  de_genes_in_module <- row.names(de_results[de_results["padj"]<=0.05,][row.names(de_results[de_results["padj"]<=0.05,])%in%module_genes[[module]],])
  num_de_genes_in_module <- length(de_genes_in_module)
  print(paste(de_genes_in_module, collapse=','))
  overview_df<-rbind(overview_df,c('module'=module,
                                   'number of module genes'=num_module_genes,
                                   'number of DE genes'= length(row.names(de_results[de_results["padj"]<=0.05,])),
                                   'number of DE genes in module'=num_de_genes_in_module,
                                   'ratio of DE genes to module genes'=num_de_genes_in_module/num_module_genes,
                                   'DE file'=ll,
                                   'DE genes in module'=paste(de_genes_in_module, collapse=','))
  )
}
colnames(overview_df)<-c('module',
                         'number of module genes',
                         'number of DE genes',
                         'number of DE genes in module',
                         'ratio of DE genes to module genes',
                         'DE file',
                         'DE genes in module')
overview_df<-column_to_rownames(overview_df, 'module')
write.csv(overview_df,file.path(file.path(params[["output_dir"]],'overview.csv')))