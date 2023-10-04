library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)
library(pathview)

make_go_enr = function(genelist, 
                       fname, 
                       folder, 
                       fromType="SYMBOL",
                       OrgDb=org.Mm.eg.db,
                       prefixes=c("BP","MF","CC","All"),
                       image_width=30,
                       image_height=30){

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

        p1 = dotplot(ego, showCategory=30) + ggtitle(paste0("GO enrichment ", prefix))
        ggsave(p1, filename = file.path(folder,fname, paste0("dotplot_go_", prefix,"_",fname,".jpg")), width = image_width, height=image_height, units="cm")
  
        p2 = barplot(ego, showCategory=30) + ggtitle(paste0("GO enrichment ", prefix))
        ggsave(p2, filename = file.path(folder, fname, paste0("barplot_go_", prefix,"_",fname,".jpg")), width = image_width, height=image_height, units="cm")
      }
  }
}

make_kegg_enr <- function(deg_input,
                          fname,
                          folder,
                          OrgDb=org.Mm.eg.db,
                          pvalueCutoff=0.05,
                          FC_col="log2FoldChange",
                          degs=NULL){
  
  old_wd <- getwd()
  dir.create(file.path(folder, 'kegg'))
    
  if(species(OrgDb)=="Mus musculus"){
        organism="mmu"
  } else if (species(OrgDb)=="Homo sapiens"){
        organism="hsa"
  } else {stop(paste0("Not a valid organism! ",species(OrgDb)))}
  
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

create_geneModulePathwayFile <- function(modules,root_filepath,prefixes=c("BP","MF","CC","All")){
    #print(modules)
    for (prefix in prefixes){
        multiple_modules<-data.frame("Gene"=c(),"Module"=c(),"Pathway"=c(),"Pathway ID"=c())
        for (module in names(modules)){
            if(file.exists(file.path(root_filepath,module,paste0("go_enr_",prefix,"_",module,".csv")))){
                go_file<-read.csv(file.path(root_filepath,module,paste0("go_enr_",prefix,"_",module,".csv")))
                #print("pathway_vec")
                #pathways<-sapply(modules[[module]],function(gene){go_file[grep(gene,go_file[,"geneID"],c("Description","ID")]})
                #print(pathways)                                                 
                #print(paste0(go_file[grep(modules[[module]][1],go_file[,"geneID"]),"ID"],collapse=','))
                pathway_vec<-sapply(modules[[module]],function(gene,go_file){paste0(go_file[grep(gene,go_file[,"geneID"]),"Description"],collapse=',')},go_file)
                #print("pathwayID_vec")
                pathwayID_vec<-sapply(modules[[module]],function(gene,go_file){paste0(go_file[grep(gene,go_file[,"geneID"]),"ID"],collapse=',')},go_file)
                single_module<-data.frame("Gene"=modules[[module]],"Module"=module,"Pathway"=pathway_vec,"Pathway ID"=pathwayID_vec)
                #print(single_module)
                #multiple_modules<-c(multiple_modules,single_module)
                multiple_modules<-rbind(multiple_modules,single_module)
                #print(multiple_modules)
                }
            }
        
        write.csv(multiple_modules,file.path(root_filepath,paste0("geneModulePathway_",prefix)))
        }
    #stop("End of function: create_geneModulePathwayFile")
}

return_literatureModuleGeneOverview<-function(
    module_entities, 
    literatureDF, 
    entity_col,
    pupId_col)
{
    #This function uses a text mining generated file containing found literature to molecular entities, 
    #such as genes or proteins and returns a data.frame that gives an overview about found literature for entities contained 
    #in WGCNA computed modules.
    #
    #module_entities should be a named list of modules containing molecular entities.
    #literatureDF is a data.frame containing literature results containg molecular entities
    #entity_col is the column identifier containing entity identifiers in literatureDF.
    #pupId_col is the column identifier containing publication identifiers in literatureDF
    
    #Check if module gene is a list with content
    stopifnot(is.list(module_entities) && !identical(module_entities, list()))
    #Check if literature DF is a dataframe -> Can be empty, whereas module genes should not be empty.
    stopifnot(is.data.frame(literatureDF))
    
    #Initiate overview data frame
    literatureModuleGeneOverview<-data.frame()
    
    for (module in names(module_entities)){
        
        intersection_selection<-literatureDF[,entity_col]%in% module_entities[[module]]
        module_literature<-literatureDF[intersection_selection,]
        
        literatureModuleGeneOverview<-rbind(literatureModuleGeneOverview,c("module"=module, 
                      "number of module genes"=length(module_entities[[module]]),
                      "number of module gene publications"=dim(module_literature)[1],
                      "found module genes in literature"=paste(module_literature[!duplicated(module_literature[,entity_col]),entity_col], collapse=','),
                      "publication IDs"=paste(module_literature[,pupId_col],collapse=',')))
        
    }
    colnames(literatureModuleGeneOverview)<-c("module", 
                      "number of module genes",
                      "number of module gene publications",
                      "found module genes in literature","publication IDs")
    
    #return results
    literatureModuleGeneOverview
}