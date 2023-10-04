library(tidyverse)
library(jsonlite)
library(org.Mm.eg.db)

library(WGCNA)
library(flashClust)

source(file.path("src","shared","WGCNA","enrichment_WGCNA_functions.R"))

perform_WGCNA <- function(datExpr,datTraits,result_folder,de_dir,literatureDF,als_markers,mapping=NULL,image_width=30,OrgDb=org.Mm.eg.db){
    
    dir.create(result_folder)
    result_folder_name<-rev(strsplit(result_folder, split ="/")[[1]])[1]
    
    datExpr = t(datExpr)  # now samples are rows and genes are columns
    
    # Run this to check if there are gene outliers
    gsg = goodSamplesGenes(datExpr, verbose = 3)
    #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
    if (!gsg$allOK){
        if (sum(!gsg$goodGenes)>0){
            printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
            if (sum(!gsg$goodSamples)>0){
                printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
                datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
            }
        }
    }
    #should return TRUE if datasets align correctly, otherwise your names are out of order
    stopifnot(all(rownames(datTraits)==rownames(datExpr)))

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
    
    jpeg(file=file.path(result_folder,"DendroAndColors.jpg"))
    plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")
    dev.off()
    #ggsave(DendroandColors, filename = "DendroandColors.jpg", width = 20, units="cm")

    # Choose a soft threshold power
    powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
    sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
    pdf(file=file.path(result_folder,"pickSoftThreshold.pdf"),width=9,height=5)
    par(mfrow= c(1,2))
    cex1=0.9
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))+
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")+
    abline(h=0.90, col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))+
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
    dev.off()
    #ggsave(sftp_sftmf, filename = "sftp_sftmf.jpg", width = 20, units="cm")
    #ggsave(sftp_mc, filename = "sftp_mc.jpg", width = 20, units="cm")
    
    #from this plot choose the power for which the scale free topology index reaches 0.90
    if (any((-sign(sft$fitIndices[,3])*sft$fitIndices[,2])>=0.9)){
        softPower = sft$fitIndices[(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])>=0.9,1][1]
        }
    else{
        #print(sft$fitIndices[,1])
        #print(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
        #print(which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2])==max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])))
        softPower = sft$fitIndices[which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2])==max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])),1]
        #print(softPower)
        #stop("debugging softPower")
        }
    
    enableWGCNAThreads()
    # see scale independence file "softthreshold_SOD1.png"
    adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
    arrow::write_feather(as.data.frame(adjacency), file.path(params[["output_dir"]],"adjacency_matrix.feather"))
    # Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
    
    #translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
    TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
    dissTOM = 1-TOM
    
    # Generate Modules --------------------------------------------------------
    
    # Generate a clustered gene tree
    geneTree = flashClust(as.dist(dissTOM), method="average")
    
    jpeg(file=file.path(result_folder,"geneTree.jpg"))
    plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
    dev.off()
    
    #This sets the minimum number of genes to cluster into a module
    minModuleSize = 30 
    dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
    
    dynamicColors= labels2colors(dynamicMods)
    MEList= moduleEigengenes(datExpr, colors= dynamicColors)#,softPower = 18)
    MEs= MEList$eigengenes
    MEDiss= 1-cor(MEs)
    METree= flashClust(as.dist(MEDiss), method= "average")

    #save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
    
    #plots tree showing how the eigengenes cluster together
    jpeg(file=file.path(result_folder,"MEtree.jpg"))
    plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
    dev.off()
    
    #set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
    MEDissThres = 0.25
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
    mergedColors = merge$colors
    mergedMEs = merge$newMEs
    
    #plot dendrogram with module colors below it
    jpeg(file=file.path(result_folder,"DendroAndColors.jpg"))
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
    dev.off()
    
    moduleColors = mergedColors
    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(moduleColors, colorOrder)-1
    MEs = mergedMEs
    
    #save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

    # Correlate traits --------------------------------------------------------
    
    #Define number of genes and samples
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    
    #Recalculate MEs with color labels
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, datTraits, use= "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    write.csv(moduleTraitCor,file.path(result_folder,paste0("moduleTraitCor.csv")))
    
    
    #Print correlation heatmap between modules and traits
    jpeg(file=file.path(result_folder,"corrHeatmapAll.jpg"))
    textMatrix= paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep= "")
    dim(textMatrix)= dim(moduleTraitCor)
    par(mar=  c(6, 3, 3, 3))
    
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

    if(species(OrgDb)=="Mus musculus"){
        condition="mut"
    } else if (species(OrgDb)=="Homo sapiens"){
        condition="als"
    } else {stop(paste0("Not a valid organism! ",species(OrgDb)))}
    
    #gets modules where the mut conditions are highly correleted
    module_selection<-(abs(moduleTraitCor[,"ctrl_female"])<0.5&abs(moduleTraitCor[,"ctrl_male"])<0.5)&(abs(moduleTraitCor[,paste0(condition,"_female")])>0.5|abs(moduleTraitCor[,paste0(condition,"_male")])>0.5)

    #check if is dataframe, otherwise it might be a single vector containing 1> modules in selection
    if(is.data.frame(moduleTraitCor[module_selection,])){
        
         #display the corelation values with a heatmap plot
         jpeg(file=file.path(result_folder,"corrHeatmapSelected.jpg"))
         par(mar=  c(6, 9, 3, 3))
         corrHeatmapSelected<-labeledHeatmap(Matrix= moduleTraitCor[module_selection,], 
                xLabels= names(datTraits), 
                yLabels= names(MEs)[module_selection], 
                ySymbols= names(MEs)[module_selection], 
                colorLabels= FALSE, 
                colors= blueWhiteRed(50), 
                textMatrix= textMatrix[module_selection,], 
                setStdMargins= FALSE, 
                cex.text= 0.5, 
                zlim= c(-1,1), 
                main= paste("Module-trait relationships"))
         dev.off()
         }
    
    #get selected genes
    all_Colors<-sapply(names(MEs),function(ME){sub("ME","",ME)})
    #selected_Colors<-sapply(names(MEs)[module_selection],function(ME){sub("ME","",ME)})
    #names(datExpr)[moduleColors=="brown"]
    module_genes<-sapply(all_Colors,function(module_color){rownames(t(datExpr)[mergedColors==module_color,])})
    saveRDS(module_genes,file=file.path(result_folder,'module_genes.rds'))
    write(toJSON(module_genes),file=file.path(result_folder,'module_genes.json'))
    
    
    ## enrichment of each module
    overview_df<- data.frame()
    
    for (module in names(module_genes)){
        print(module)
        ##GO enrichment
        make_go_enr(module_genes[[module]],module,file.path(result_folder),OrgDb=OrgDb)
        
        ## kegg enrichment
        de_symbol_col<-"gene_name"
        types = c("only_females", "only_males")#[c(moduleTraitCor[module,"mut_female"]>0.5,moduleTraitCor[module,"mut_male"]>0.5)]
        for(ll in types){
            obj <- file.path(de_dir,ll,"/mut_vs_ctrl.csv")
            de_results <- read.csv(obj)
            
            if(!is.null(mapping)){
                de_results<-merge(de_results,mapping,by="UniProtName")
                      }
            
            de_results<-na.omit(de_results[!duplicated(de_results[,de_symbol_col]), ])
            rownames(de_results) = de_results[,de_symbol_col]
            
            make_kegg_enr(de_results[rownames(de_results) %in% module_genes[[module]],],
                  paste0(module,'_',ll),
                  file.path(result_folder,module),
                  degs=na.omit(de_results[de_results["padj"]<=0.05,]),
                  OrgDb=OrgDb)
        
        ## get module overview table
        num_module_genes <- length(module_genes[[module]])
        de_genes_in_module <- intersect(de_results[de_results["padj"]<=0.05,'gene_name'],module_genes[[module]]) 
        num_de_genes_in_module <- length(de_genes_in_module)
        als_markers_in_module <- intersect(sapply(als_markers, tolower), sapply(module_genes[[module]], tolower))
            
        #print(als_markers_in_module)
        #print(sapply(als_markers, tolower))
        #print(module_genes[[module]])
        #print(sapply(module_genes[[module]], tolower))
        #print(intersect(sapply(als_markers, tolower), sapply(module_genes[[module]], tolower)))
        #stop()
    
        overview_df<-rbind(overview_df,c('module'=module, 
                      'number of module genes'=num_module_genes,
                      'number of DE genes'= length(row.names(de_results[de_results["padj"]<=0.05,])),
                      'number of DE genes in module'=num_de_genes_in_module,
                      'ratio of DE to module genes'=num_de_genes_in_module/num_module_genes,
                      'DE file'=ll,
                      'DE genes in module'=paste(de_genes_in_module, collapse=','),
                      'ALS markers in module'=paste(als_markers_in_module, collapse=','))
                          )
        }
    }
    colnames(overview_df)<-c('module',
                         'number of module genes',
                         'number of DE genes',
                         'number of DE genes in module',
                         'ratio of DE to module genes',
                         'DE file',
                         'DE genes in module',
                         'ALS markers in module')
    #Creating overview files
    write.csv(overview_df,file.path(result_folder,paste0('overview_',rev(strsplit(result_folder, split ="/")[[1]])[1],'.csv')))
    create_geneModulePathwayFile(module_genes,result_folder) #creating overview file that shows found pathways per gene
    write.csv(paste0("softPowerThreshold: ",softPower),file.path(result_folder,paste0("parameters_",result_folder_name,".csv")))
    #Creating module gene to literature overview
    write.csv(return_literatureModuleGeneOverview(module_genes, literatureDF, "term", "b_pubmed_id"),
             file.path(result_folder,paste0("literatureOverview_",result_folder_name,".csv")))
}