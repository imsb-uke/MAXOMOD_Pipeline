#!/usr/bin/env Rscript

library(MOFA2)
library(tidyverse)
library(DESeq2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(argparse)
#pathway enrichment
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
#library(MOFAdata)
#library(msigdb)
library(msigdbr)
library(reshape2)
library(pheatmap)

variance_threshold <- function(data,threshold){
    if (nrow(data)<=threshold){ret}
    variances <- apply(data,1,var)
    reduced_data <- data[variances<threshold,]
    print(paste0("Data reduced from ",dim(data)[1]," to ",dim(reduced_data)[1], " features by a variance threshold of ",threshold,"."))
    reduced_data
}

topN_byVariance <- function(data,n){
    stopifnot(n>0)
    if (nrow(data)<=n){
        data
    } else {
        variances <- apply(data,1,var)
        top_features <- names(sort(variances,decreasing = TRUE))[1:n]
        data_reduced <- data[top_features,]
        data_reduced
    }
}

get_args <- function() {

    # create parser object
    parser <- ArgumentParser()
    parser$add_argument("--dataset-name", type = "character", default = NULL, help = "dataset key from params.yaml to use")
    parser$add_argument("--mofa-type", type = "character", default = NULL, help = "mofa stage name dependent on which kind of omics where used for the analysis")

    args <- parser$parse_args()
    args
}

# MAIN

## Read params

args = get_args()

dataset_name = args$dataset_name[1]
print(dataset_name)

mofa_type = args$mofa_type[1]
print(mofa_type)

params <- yaml::read_yaml("params.yaml")[[mofa_type]][["stages"]][[dataset_name]]

## Other Parameters
TRANSCRIPTOMICS_PARAMS_NAME="transcriptomics"
PROTEOMICS_PARAMS_NAME="proteomics"
SRNA_MATURE_PARAMS_NAME="srna_mature"
SRNA_HAIRPIN_PARAMS_NAME="srna_hairpin"
METABOLOMICS_PARAMS_NAME="metabolomics"
PHOSPHOPROTEOMICS_PARAMS_NAME="phosphoproteomics"
MAX_SEED=10000
GROUP_NAME="group1"
MAX_FACTOR_INTERVAL=10

###create main directory
out_dir <- params[["out_dir"]]
dir.create(file.path(params[["out_dir"]]), recursive = TRUE)

### data import
#A list of matrices, where each entry corresponds to one view/omics. Samples are stored in columns and features in rows.

#import rnaseq and metadata
if (TRANSCRIPTOMICS_PARAMS_NAME %in% names(params[["omics"]])){
    print("rnaseq")
    rnaseq_meta <- read.csv(file.path(params[["omics"]][["transcriptomics"]][["annotation_file"]]),row.names=1)
    rnaseq_meta <- tibble::rownames_to_column(rnaseq_meta, "SampleID")
   # print(rnaseq_meta)

    gene_expression <- as.data.frame(read.csv(file.path(params[["omics"]][["transcriptomics"]][["expression_file"]]),row.names=1))
    gene_expression <- tibble::rownames_to_column(gene_expression, "gene_id") %>%
        mutate(gene_id=sub("[.][0-9]*", "", gene_id)) #get rid of version number in gene id
    #print(head(gene_expression))

    #get gene symbol from ensembl id
    if (params[["organism"]]=="Mus musculus"){
        OrgDb <- org.Mm.eg.db
        } else if (params[["organism"]]=="Homo sapiens"){
        OrgDb <- org.Hs.eg.db
        } else {
        stop(paste0(params[["organism"]]," is not an valid organism. Should be one of 'Mus musculus' or'Homo sapiens'"))
    }
    mapping <- mapIds(OrgDb, keys = pull(gene_expression,"gene_id"), column = c("SYMBOL"), keytype = c("ENSEMBL"), multiVals="first")
    mapping <- data.frame(mapping) %>%
        tibble::rownames_to_column("gene_id")
    #print(head(mapping))
    gene_expression<-merge(gene_expression, mapping, by='gene_id') %>%
        drop_na() %>%
        dplyr::rename(gene_name = mapping)

    if (params[["organism"]]=="Homo sapiens"){
        gene_expression<-gene_expression %>%
            rename_with(~ gsub("[.]", "-", .x)) #replace dot in sample names
    }

    #########get rid of duplicates, set gene name as rownames and deselect others#########
    gene_expression<-gene_expression[!duplicated(gene_expression$gene_name),] %>%
        remove_rownames %>%
        tibble::column_to_rownames(var="gene_name") %>%
        dplyr::select(-gene_id)

    #top features by variance
    #gene_expression<-topN_byVariance(gene_expression,3000)

    #print(head(gene_expression))

    gene_expression_norm <- as.data.frame(assay(vst(DESeqDataSetFromMatrix(countData = as.matrix(gene_expression),colData = rnaseq_meta,design = ~ Condition))))
    head(gene_expression_norm)

    pdf(file=file.path(out_dir,paste0(TRANSCRIPTOMICS_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(gene_expression_norm))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        gene_expression_norm <- topN_byVariance(gene_expression_norm,params[["topN_byVar"]])
    }
}

#import proteomics
if (PROTEOMICS_PARAMS_NAME %in% names(params[["omics"]])){
    print("proteomics")
    proteomics_meta <- read.csv(file.path(params[["omics"]][["proteomics"]][["annotation_file"]]),row.names=1)
    proteomics_meta <- tibble::rownames_to_column(proteomics_meta, "SampleID")
    #print(proteomics_meta)

    protein_abundances <- as.data.frame(t(read.csv(file.path(params[["omics"]][["proteomics"]][["expression_file"]]),row.names=1))%>%
        apply(2,function(x){x - mean(x)})) #mean normalize protein data

    #print(head(protein_abundances))

    pdf(file=file.path(out_dir,paste0(PROTEOMICS_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(protein_abundances))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        protein_abundances <- topN_byVariance(protein_abundances,params[["topN_byVar"]])
    }
}

#import small rna (mature)
if (SRNA_MATURE_PARAMS_NAME %in% names(params[["omics"]])){
    print("small rna mature")
    srna_mature_meta <- read.csv(file = file.path(params[["omics"]][["srna_mature"]][["annotation_file"]]))
    #print(srna_mature_meta)

    mature_to_keep <- readRDS(file.path(params[["omics"]][["srna_mature"]][["keep_file"]]))
    srna_mature_expression <- as.data.frame(readRDS(file = file.path(params[["omics"]][["srna_mature"]][["expression_file"]]))[mature_to_keep,])
    #print(head(srna_mature_expression))

    srna_mature_norm <- as.data.frame(assay(varianceStabilizingTransformation(DESeqDataSetFromMatrix(countData = as.matrix(srna_mature_expression),
                           colData = srna_mature_meta,
                           design = ~ Condition))))
    head(srna_mature_norm)

    pdf(file=file.path(out_dir,paste0(SRNA_MATURE_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(srna_mature_norm))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        srna_mature_norm <- topN_byVariance(srna_mature_norm,params[["topN_byVar"]])
    }
}

#import small rna (hairpin)
if (SRNA_HAIRPIN_PARAMS_NAME %in% names(params[["omics"]])){
    print("small rna hairpin")
    srna_hairpin_meta <- read.csv(file = file.path(params[["omics"]][["srna_hairpin"]][["annotation_file"]]))
    #print(srna_hairpin_meta)

    hairpin_to_keep <- readRDS(file.path(params[["omics"]][["srna_hairpin"]][["keep_file"]]))
    srna_hairpin_expression <- as.data.frame(readRDS(file = file.path(params[["omics"]][["srna_hairpin"]][["expression_file"]]))[hairpin_to_keep,])
    #print(head(srna_hairpin_expression))

    srna_hairpin_norm <-  as.data.frame(assay(varianceStabilizingTransformation(DESeqDataSetFromMatrix(countData = as.matrix(srna_hairpin_expression),
                           colData = srna_hairpin_meta,
                           design = ~ Condition))))
    head(srna_hairpin_norm)

    pdf(file=file.path(out_dir,paste0(SRNA_HAIRPIN_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(srna_hairpin_norm))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        srna_hairpin_norm <- topN_byVariance(srna_hairpin_norm,params[["topN_byVar"]])
    }
}

#import metabolomics
if (METABOLOMICS_PARAMS_NAME %in% names(params[["omics"]])){
    print("metabolomics")
    print("Importing files")
    metabolomics_meta <- read.csv(file = file.path(params[["omics"]][[METABOLOMICS_PARAMS_NAME]][["annotation_file"]]))
    #print(metabolomics_meta)

    metabolomics_abundance <- as.data.frame(read.csv(file = file.path(params[["omics"]][[METABOLOMICS_PARAMS_NAME]][["expression_file"]]),row.names=1))

    if (!is.null(params[["omics"]][[METABOLOMICS_PARAMS_NAME]][["sample_mapping"]])){
        metabolomics_mapping <- read.delim(file = file.path(params[["omics"]][[METABOLOMICS_PARAMS_NAME]][["sample_mapping"]]))
        #print(metabolomics_mapping)

        metabolomics_aundance <- metabolomics_abundance[,metabolomics_mapping$MetSampleID]
        names(metabolomics_abundance) <- metabolomics_mapping$SampleID
    }

    #print(head(metabolomics_abundance))

    pdf(file=file.path(out_dir,paste0(METABOLOMICS_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(metabolomics_abundance))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        metabolomics_abundance <- topN_byVariance(metabolomics_abundance,params[["topN_byVar"]])
    }
}
#import phosphoproteomics
if (PHOSPHOPROTEOMICS_PARAMS_NAME %in% names(params[["omics"]])){
    print("phosphoproteomics")
    phosphoproteomics_meta <- read.csv(file.path(params[["omics"]][["phosphoproteomics"]][["annotation_file"]]),row.names=1)
    phosphoproteomics_meta <- tibble::rownames_to_column(phosphoproteomics_meta, "SampleID")
    #print(phosphoproteomics_meta)

    phosphoproteomics_abundance <- as.data.frame(t(read.csv(file.path(params[["omics"]][["phosphoproteomics"]][["expression_file"]]),row.names=1)))
    #print(head(phosphoproteomics_abundance))
    ### Change sample names by mapping -> to include phospho data from separate samples -> do not use in general -> in this case the argument is that the mouse data is so similar that it won't have too much of an effect, even if the data comes from different mice ###
    if (!is.null(params[["omics"]][["phosphoproteomics"]][["sample_matching"]])){
        phospho_sample_matching <- readr::read_tsv(params[["omics"]][["phosphoproteomics"]][["sample_matching"]])

        phospho_sample_matching<-phospho_sample_matching %>%
            column_to_rownames("AdditionalSampleID")
        #print(phospho_sample_matching)

        colnames(phosphoproteomics_abundance) <- phospho_sample_matching[colnames(phosphoproteomics_abundance),c("OriginalSampleID")]
        #print(head(phosphoproteomics_abundance))

        phosphoproteomics_meta$SampleID<-phospho_sample_matching[phosphoproteomics_meta$SampleID,c("OriginalSampleID")]
        #print(phosphoproteomics_meta)
    }

    pdf(file=file.path(out_dir,paste0(PHOSPHOPROTEOMICS_PARAMS_NAME,"_variancePerFeature.pdf")))
    var_hist <- hist(log(rowVars(as.matrix(phosphoproteomics_abundance))),breaks=100)
    var_hist$counts <- cumsum(var_hist$counts)    # Change histogram counts
    plot(var_hist)
    dev.off()

    if (!is.null(params[["topN_byVar"]])) {
        phosphoproteomics_abundance <- topN_byVariance(phosphoproteomics_abundance,params[["topN_byVar"]])
    }

}

#Collect omics data and meta data into lists
omics_data <- list()
omics_meta <- list()
for (omic_type in names(params[["omics"]])){
    if (omic_type==TRANSCRIPTOMICS_PARAMS_NAME){
        omics_data[[TRANSCRIPTOMICS_PARAMS_NAME]]<-gene_expression_norm
        omics_meta[[TRANSCRIPTOMICS_PARAMS_NAME]]<-rnaseq_meta
    }
    if (omic_type==PROTEOMICS_PARAMS_NAME){
        omics_data[[PROTEOMICS_PARAMS_NAME]]<-protein_abundances
        omics_meta[[PROTEOMICS_PARAMS_NAME]]<-proteomics_meta
    }
    if (omic_type==SRNA_MATURE_PARAMS_NAME){
        omics_data[[SRNA_MATURE_PARAMS_NAME]]<-srna_mature_norm
        omics_meta[[SRNA_MATURE_PARAMS_NAME]]<-srna_mature_meta
    }
    if (omic_type==SRNA_HAIRPIN_PARAMS_NAME){
        omics_data[[SRNA_HAIRPIN_PARAMS_NAME]]<-srna_hairpin_norm
        omics_meta[[SRNA_HAIRPIN_PARAMS_NAME]]<-srna_hairpin_meta
    }
    if (omic_type==METABOLOMICS_PARAMS_NAME){
        omics_data[[METABOLOMICS_PARAMS_NAME]]<-metabolomics_abundance
        omics_meta[[METABOLOMICS_PARAMS_NAME]]<-metabolomics_meta
    }
    if (omic_type==PHOSPHOPROTEOMICS_PARAMS_NAME){
        omics_data[[PHOSPHOPROTEOMICS_PARAMS_NAME]]<-phosphoproteomics_abundance
        omics_meta[[PHOSPHOPROTEOMICS_PARAMS_NAME]]<-phosphoproteomics_meta
    }
}
#print(omics_meta)

#Check dimensions
print("Omic data dimensions:")
for(omic_data in omics_data){
    print(dim(as.matrix(omic_data)))
    }

### MAIN MOFA execution ###

for(sex_selection in params[["sex_selection"]]){
    #Find shared column names
    #shared_cols <- Reduce(intersect, lapply(omics_data,colnames))

    #print("Shared samples in omics data:")
    #print(shared_cols)

    #join all meta data
    shared_meta_cols <- Reduce(intersect, lapply(omics_meta,colnames))

    print("Shared columns in omics meta data:")
    print(shared_meta_cols)

    joined_meta <- Reduce(function(omics_meta1,omics_meta2){full_join(omics_meta1,omics_meta2,by=shared_meta_cols)}, omics_meta)
    #joined_meta <- joined_meta %>%
    #    filter(SampleID %in% shared_cols)
    #print(joined_meta)

    #Checks if sample names (colnames) from omic data in joined_meta and then arranges the cols after joined_meta
    for(omic_data in omics_data){
        stopifnot(all(colnames(omic_data %in% joined_meta$SampleID)))
        omic_data <- omic_data[,joined_meta$SampleID[joined_meta$SampleID %in% colnames(omic_data)]]
    }

    #model metadata needs a column called 'sample'
    joined_meta <- joined_meta %>%
        dplyr::rename(sample=SampleID)

    #print("joined")
    #print(joined_meta)

    if(sex_selection!="all"){
        out_dir <- file.path(params[["out_dir"]],sex_selection)
        dir.create(out_dir, recursive = TRUE)

        joined_meta <- joined_meta %>%
            dplyr::filter(Sex==sex_selection)
        #print(joined_meta)
    }

    MOFAobject <- create_mofa(lapply(omics_data,function(omic_data){
        for(missingSample in dplyr::setdiff(joined_meta$sample,colnames(omic_data))){
            omic_data[missingSample] <- rep(c(NA),times=nrow(omic_data))}
        as.matrix(omic_data[,joined_meta$sample])}))

    print(MOFAobject)

    ggsave(file.path(out_dir,"data_overview.pdf"),
    plot_data_overview(MOFAobject)
    )

    #Define data options
    #scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE
    #scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE

    data_opts <- get_default_data_options(MOFAobject)
    data_opts$scale_views<-TRUE
    head(data_opts)

    #Define model options
    #num_factors: number of factors
    #likelihoods: likelihood per view (options are “gaussian”, “poisson”, “bernoulli”). By default they are learnt automatically. We advise users to use “gaussian” whenever possible!
    #spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE.
    #spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE.
    #ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
    #ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.

    #Only change the default model options if you are familiar with the underlying mathematical model!

    model_opts <- get_default_model_options(MOFAobject)

    head(model_opts)

    #Prepare the MOFA object
    if (is.null(params[["model_training"]][["seed"]])){
        if (!is.null(params[["model_training"]][["n_initializations"]])) {
            n_inits <- params[["model_training"]][["n_initializations"]]
        } else {
            n_inits <- 10
        }

        mofa_models <- list()
        for (i in 1:n_inits){
            #Define train options
            #maxiter: number of iterations. Default is 1000.
            #convergence_mode: “fast”, “medium”, “slow”. For exploration, the fast mode is good enough.
            #startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence).
            #freqELBO: frequency of computations of the ELBO.
            #gpu_mode: use GPU mode? (needs cupy installed and a functional GPU).
            #stochastic: use stochastic inference? (default is FALSE).
            #verbose: verbose mode?
            #seed: random seed

            train_opts <- get_default_training_options(MOFAobject)
            train_opts$seed <- sample(1:MAX_SEED,1)
            #train_opts$convergence_mode <- "slow"
            #train_opts$maxiter <- 5000
            head(train_opts)

            MOFAobject <- prepare_mofa(
              object = MOFAobject,
              data_options = data_opts,
              model_options = model_opts,
              training_options = train_opts
            )

            #Train the MOFA model

            outfile = file.path(out_dir,paste0("mofa_model_seed",train_opts$seed,".hdf5"))
            MOFAobject.trained <- run_mofa(MOFAobject, outfile)
            file.rename(from=outfile,
                        to=file.path(out_dir,paste0("mofa_model_seed",train_opts$seed,
                                                    "_ELBO", get_elbo(MOFAobject.trained), ".hdf5")))
            print(train_opts$seed)
            mofa_models[[toString(train_opts$seed)]] <- MOFAobject.trained
            print(typeof(MOFAobject.trained))
        }
        print(mofa_models)
        MOFAobject.trained <- select_model(mofa_models)
        compare_elbo(mofa_models, log = TRUE)
        ggsave(file.path(out_dir,"compare_elbo.pdf"))
        pdf(file=file.path(out_dir,"compare_factors.pdf"))
        compare_factors(mofa_models)
        dev.off()

    } else {
        #Define train options
        #maxiter: number of iterations. Default is 1000.
        #convergence_mode: “fast”, “medium”, “slow”. For exploration, the fast mode is good enough.
        #startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence).
        #freqELBO: frequency of computations of the ELBO.
        #gpu_mode: use GPU mode? (needs cupy installed and a functional GPU).
        #stochastic: use stochastic inference? (default is FALSE).
        #verbose: verbose mode?
        #seed: random seed

        train_opts <- get_default_training_options(MOFAobject)
        train_opts$seed <- params[["model_training"]][["seed"]]
        #train_opts$drop_factor_threshold <- 0.02 #Uncomment to drop factors with less than 2% explained variance, Otherwise mofa is initiated with 15 factors and factors with explained variance of 0 are dropped
        #train_opts$convergence_mode <- "slow"
        #train_opts$maxiter <- 5000
        head(train_opts)

        MOFAobject <- prepare_mofa(
              object = MOFAobject,
              data_options = data_opts,
              model_options = model_opts,
              training_options = train_opts
        )

        #Train the MOFA model

        outfile = file.path(out_dir,paste0("mofa_model_seed",params[["model_training"]][["seed"]],".hdf5"))
        MOFAobject.trained <- run_mofa(MOFAobject, outfile)
        file.rename(from=outfile,
                    to=file.path(out_dir,
                                 paste0("mofa_model_seed",params[["model_training"]][["seed"]],
                                        "_ELBO", get_elbo(MOFAobject.trained), ".hdf5")))
    }

    seed <- get_default_training_options(MOFAobject.trained)$seed

    #Add metadata to the model

    #The metadata is stored as a data.frame object in model@samples_metadata, and it requires at least the column sample. The column group is required only if you are doing multi-group inference.
    #The number of rows must match the total number of samples in the model (sum(model@dimensions$N)).

    samples_metadata(MOFAobject.trained) <- joined_meta
    head(MOFAobject.trained@samples_metadata, n=5)

    pdf(file=file.path(out_dir,paste0(seed,"_factor_correlation.pdf"))
    )
    plot_factor_cor(MOFAobject.trained)
    dev.off()


    #Variance decomposition
    print("Variance decomposition")
    variance_explained <- get_variance_explained(MOFAobject.trained)
    print(variance_explained)

    for (expl_var_name in names(variance_explained)){
        write.csv(variance_explained[[expl_var_name]][[GROUP_NAME]],
                  file.path(out_dir,paste0(seed,"_",expl_var_name,".csv")))
    }

    plot_variance_explained(MOFAobject.trained, x="view", y="factor")
    ggsave(file.path(out_dir,paste0(seed,"_variance_decomposition.pdf")))

    plot_variance_explained(MOFAobject.trained, x="view", y="factor", plot_total = T)
    ggsave(file.path(out_dir,paste0(seed,"_variance_decomposition_total.pdf")))

    # Fetch factors and save as csv
    factor_matrix <- get_factors(MOFAobject.trained, as.data.frame = FALSE)[[GROUP_NAME]]
    write.csv(factor_matrix,
              file.path(out_dir,paste0(seed,"_factor_matrix.csv")))
    #print(factor_matrix)

    factors <- colnames(factor_matrix)

    if (sex_selection=="all") {
        covariates <- colnames(joined_meta %>%
            dplyr::select(-sample,-CellLine))
    } else {
        covariates <- colnames(joined_meta %>%
            dplyr::select(-sample,-CellLine,-Sex))
    }

    if((length(factors)%%MAX_FACTOR_INTERVAL)==0){
        batches<-(length(factors)/MAX_FACTOR_INTERVAL)
    }else{
        batches<-(length(factors)/MAX_FACTOR_INTERVAL)+1
    }

    for (annotation in covariates){
        print(annotation)
        for(x in 1:batches)
        {
            if ((x*MAX_FACTOR_INTERVAL)>length(factors)){
                start<-((x-1)*MAX_FACTOR_INTERVAL)+1
                end<-length(factors)
            } else {
                start<- ((x-1)*MAX_FACTOR_INTERVAL)+1
                end<-x*MAX_FACTOR_INTERVAL
            }

            if (annotation=="Sex"){
                p <- plot_factors(
                MOFAobject.trained,
                factors = start:end,
                color_by = annotation
                )+ # The output of plot_factor is a ggplot2 object that we can edit
                scale_color_manual(values=c("male"="black", "female"="red")) +
                scale_fill_manual(values=c("male"="black", "female"="red"))
            }else{
                p <- plot_factors(
                MOFAobject.trained,
                factors = start:end,
                color_by = annotation
            )}
            ggsave(file.path(out_dir,paste0(seed,"_factors_",annotation,"_",x,".pdf")),p)

            if (annotation=="Sex"){
                plot_factor(MOFAobject.trained,
                factors = start:end,
                color_by = annotation,
                dot_size = 3,        # change dot size
                dodge = T,           # dodge points with different colors
                legend = T,          # remove legend
                add_violin = T,      # add violin plots,
                violin_alpha = 0.25  # transparency of violin plots
                )+ # The output of plot_factor is a ggplot2 object that we can edit
                scale_color_manual(values=c("male"="black", "female"="red")) +
                scale_fill_manual(values=c("male"="black", "female"="red"))
            }else{
                plot_factor(MOFAobject.trained,
                factors = start:end,
                color_by = annotation,
                dot_size = 3,        # change dot size
                dodge = T,           # dodge points with different colors
                legend = T,          # remove legend
                add_violin = T,      # add violin plots,
                violin_alpha = 0.25  # transparency of violin plots
                )
            }
            ggsave(file.path(out_dir,paste0(seed,"_single_factors_",annotation,"_",x,".pdf")))
        }
    }

    print(params[["factors"]])
    if(!is.null(params[["factor_combinations"]])){
        for (factor_combination in params[["factor_combinations"]]){
            plot_factors(MOFAobject.trained,
                factors = factor_combination,
                color_by = "Condition",
                shape_by = "Sex"
                )
            ggsave(file.path(out_dir,paste0(seed,"_factor_combination_",paste(factor_combination,collapse='_'),".pdf")))
        }
    }

try(
ggsave(file.path(out_dir,paste0(seed,"_covariate_correlation.pdf")),
    correlate_factors_with_covariates(MOFAobject.trained,
      covariates = covariates,
      plot="log_pval"
    ))
)

#Visualisation of feature weights
#The weights provide a score for how strong each feature relates to each factor. Features with no association with the factor have values close to zero, while features with strong association with the factor have large absolute values. The sign of the weight indicates the direction of the effect: a positive weight indicates that the feature has higher levels in the cells with positive factor values, and vice versa.

#change joined_meta to heatmap annotation
annotation_samples <- tibble::column_to_rownames(joined_meta, "sample")

if (dataset_name=="human"){
    annotation_samples<-annotation_samples %>%
    mutate(age_group=cut(AgeAtDeath, breaks = c(-Inf,50,70,Inf),labels=c("<50","50-70",">70")))

    #print(annotation_samples)
}

p#rint(annotation_samples)




for (omics in views_names(MOFAobject.trained)){

    # Fetch feature matrix and save as csv
    feature_matrix <- as.data.frame(get_weights(MOFAobject.trained, view=omics, as.data.frame = FALSE)) %>%
        rename_with(~ str_replace(.x, paste0(omics,'.'),''))
    write.csv(feature_matrix,
              file.path(out_dir,paste0(seed,"_feature_matrix_",omics,".csv")))
    #print(head(feature_matrix))

    for (factor in 1:length(factors)){
        print(paste0(omics,'_',factor))
        ggsave(
            file.path(out_dir,paste0(seed,"_weights_",factor,"_",omics,".pdf")),
            plot_weights(MOFAobject.trained,
              view = omics,
              factor = factor,
              nfeatures = 10,     # Number of features to highlight
              scale = T,          # Scale weights from -1 to 1
              abs = F             # Take the absolute value?
            )
        )
        ggsave(
            file.path(out_dir,paste0(seed,"_topWeights_",factor,"_",omics,".pdf")),
            plot_top_weights(MOFAobject.trained,
              view = omics,
              factor = factor,
              nfeatures = 10
            )
        )

        #Samples sorted according to the factor values by plot_data_heatmap (decreasing)
        ggsave(
        file.path(out_dir,paste0(seed,"_heatmap_",factor,"_",omics,".pdf")),
        plot_data_heatmap(MOFAobject.trained,
              view = omics,
              factor = factor,
              features = 20,          # number of features to plot (they are selected by weight)
              annotation_samples=annotation_samples,
              # extra arguments that are passed to the `pheatmap` function
              cluster_rows = TRUE, cluster_cols = FALSE,
              show_rownames = TRUE, show_colnames = TRUE
            )
        )

        #samples clustered by pheatmap
        ggsave(
        file.path(out_dir,paste0(seed,"_heatmap_",factor,"_",omics,"_clusteredCols.pdf")),
        plot_data_heatmap(MOFAobject.trained,
              view = omics,
              factor = factor,
              features = 20,          # number of features to plot (they are selected by weight)
              annotation_samples=annotation_samples,
              # extra arguments that are passed to the `pheatmap` function
              cluster_rows = TRUE, cluster_cols = TRUE,
              show_rownames = TRUE, show_colnames = TRUE
            )
        )

        ordered_weights <- dplyr::arrange(get_weights(MOFAobject.trained,view=omics,factor=factor,as.data.frame=TRUE),desc(abs(value)))
        write.csv(ordered_weights,
                  file.path(out_dir,paste0(seed,"_orderedWeights_",factor,"_",omics,".csv")),
                  row.names=FALSE)
        tmp_topfeatures<-str_replace(as.character(droplevels(ordered_weights[1:20,c("feature")])),paste0('_(',paste(names(omics_data),collapse='|'),')'),"") #get top 20 features by absolute value and remove



        #Creates heatmap where samples ordered by condition
        tmp_names_used = intersect(rownames(annotation_samples), colnames(omics_data[[omics]]))
        ggsave(
            file.path(out_dir,paste0(seed,"_heatmap_",factor,"_",omics,"_orderedCondition.pdf")),
            pheatmap(omics_data[[omics]][tmp_topfeatures,][,tmp_names_used],
                     annotation_col=annotation_samples[tmp_names_used,], #only samples names needed
                     cluster_cols = F)
        )

        for(annotation in covariates){
            ggsave(
            file.path(out_dir,paste0(seed,"_scatter_",factor,"_",omics,"_",annotation,".pdf")),
            plot_data_scatter(MOFAobject.trained,
                  view = omics,
                  factor = factor,           # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = annotation
                )
            )
        }
    }
}

#Run umap
set.seed(seed)
DIMRED_DOTSIZE <- 3
umap <- run_umap(MOFAobject.trained,
                n_neighbors=5)
#Plot non-linear dimensionality reduction
write.csv(umap@dim_red$UMAP,
          file.path(out_dir,paste0(seed,"_umap.csv")),
          row.names=FALSE)
for(annotation in covariates){
    ggsave(file.path(out_dir,paste0(seed,"_umap_",annotation,".pdf")),
            plot_dimred(umap,
              method = "UMAP",  # method can be either "TSNE" or "UMAP"
              color_by = annotation,
              dot_size = DIMRED_DOTSIZE
            )
    )
}

    #replace protein string _MOUSE or _HUMAN
if (params[["organism"]]=="Mus musculus"){
            h_gene_sets = msigdbr(species = "mouse", category = "C2",subcategory = "CP:KEGG")
            msigdbr_df = data.frame( gene_symbol = h_gene_sets$human_gene_symbol, pathway = h_gene_sets$gs_name , dummy = rep(1, dim(h_gene_sets)[1]))

    } else if (params[["organism"]]=="Homo sapiens"){
            h_gene_sets = msigdbr(species = "human", category = "C2",subcategory = "CP:KEGG")
            msigdbr_df = data.frame( gene_symbol = h_gene_sets$human_gene_symbol, pathway = h_gene_sets$gs_name , dummy = rep(1, dim(h_gene_sets)[1]))
    } else {
    stop(paste0(params[["organism"]]," is not an valid organism. Should be one of 'Mus musculus' or'Homo sapiens'"))
}
    h_gene_sets = msigdbr(species = "human", category = "C2",subcategory = "CP:KEGG")
    msigdbr_df = data.frame( gene_symbol = h_gene_sets$human_gene_symbol, pathway = h_gene_sets$gs_name , dummy = rep(1, dim(h_gene_sets)[1]))


    msigdbr_df_mat = dcast(msigdbr_df,  pathway ~ gene_symbol, value.var="dummy", fun.aggregate=sum)
    rownames(msigdbr_df_mat) = msigdbr_df_mat$pathway
    msigdbr_df_mat = msigdbr_df_mat[,-1]

    msigdbr_df_mat[ msigdbr_df_mat > 1 ] = 1
    head(msigdbr_df_mat)

    as.matrix(msigdbr_df_mat)
    features_names(MOFAobject.trained)[[TRANSCRIPTOMICS_PARAMS_NAME]] <- toupper(features_names(MOFAobject.trained)[[TRANSCRIPTOMICS_PARAMS_NAME]])
    head(features_names(MOFAobject.trained)[[TRANSCRIPTOMICS_PARAMS_NAME]])

    #replace protein string _MOUSE or _HUMAN
if (params[["organism"]]=="Mus musculus"){
        features_names(MOFAobject.trained)[[PROTEOMICS_PARAMS_NAME]] <- sapply(features_names(MOFAobject.trained)[[PROTEOMICS_PARAMS_NAME]],str_replace,pattern="_MOUSE", replacement="")
    } else if (params[["organism"]]=="Homo sapiens"){
        features_names(MOFAobject.trained)[[PROTEOMICS_PARAMS_NAME]] <- sapply(features_names(MOFAobject.trained)[[PROTEOMICS_PARAMS_NAME]],str_replace,pattern="_HUMAN", replacement="")
    } else {
    stop(paste0(params[["organism"]]," is not an valid organism. Should be one of 'Mus musculus' or'Homo sapiens'"))
}
    head(features_names(MOFAobject.trained)[[PROTEOMICS_PARAMS_NAME]])

    for (omics in c(TRANSCRIPTOMICS_PARAMS_NAME)){#, PROTEOMICS_PARAMS_NAME)){
    enrichment.parametric <- run_enrichment(MOFAobject.trained,
                                        view = omics, factors = "all",
                                        feature.sets = as.matrix(msigdbr_df_mat),
                                        sign = "all",
                                        statistical.test = "parametric"
    )
    #print(enrichment.parametric)

    write.csv(enrichment.parametric$pval.adj,file.path(out_dir,paste0(seed,"_mofa_enrichment_",omics,"_",sex_selection,"_pvalAdj.csv")))
#    for (sigPathway in enrichment.parametric$sigPathways){
#write.csv(enrichment.parametric$pval.adj[[sigPathway]],file.path(out_dir,paste0("mofa_enrichment_rna_",sex_selection,"_",sigPathway,".csv")))
#    }

    for(factor in 1:length(factors)){

                    #write.csv(enrichment.parametric$pval.adj[[factor]],file.path(out_dir,paste0(seed,"_mofa_enrichment_",omics,"_",sex_selection,"_",factor,".csv")))

        pdf(file.path(out_dir,paste0(seed,"_mofa_enrichment_",omics,"_",sex_selection,"_",factor,".pdf")), height = 15)
        print(
        try(plot_enrichment(enrichment.parametric,
                factor = factor,
                max.pathways = 20
        )
        )
        )

        dev.off()

        pdf(file.path(out_dir,paste0(seed,"_mofa_enrichment_detailed_",omics,"_",sex_selection,"_",factor,".pdf")), height = 10)
        print(
        try(plot_enrichment_detailed(enrichment.parametric,
                         factor = factor,
                         max.genes = 8,
                         max.pathways = 10
        )
        )
        )

        dev.off()
    }
}
}