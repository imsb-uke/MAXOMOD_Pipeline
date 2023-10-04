#!/usr/bin/env Rscript

library(MOFA2)
library(argparse)
library(tidyverse)

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

params <- yaml::read_yaml("params.yaml")[["mofa_featurelist"]][["stages"]][[dataset_name]]
dir.create(file.path(params[["out_dir"]]), recursive = TRUE) 

print(params[["model"]])
# Using an existing trained model on simulated data
#file <- system.file(params[["model"]], package = "MOFA2")
model <- load_model(file.path(params[["model"]]))

# Fetch weights in data.frame format
weights <- get_weights(model, views=params[["views"]], factors=params[["factors"]], as.data.frame = TRUE)

print(head(weights))

top_features <- weights %>%
    group_by(factor,view) %>%
    #arrange(desc(abs(value))) %>%
    top_n(params[["top_n"]],wt=abs(value))

print(top_features)

group_keys <- group_keys(top_features)
print(group_keys(top_features))

groups <- group_split(top_features)

for (group_id in 1:dim(group_keys)[1]) {
    factor <- group_keys$factor[group_id]
    view <- group_keys$view[group_id]
    group <- groups[[group_id]] %>%
        arrange(desc(abs(value)))
    
    print(factor)
    print(view)
    print(group)
    
    write.table(group,
              file.path(params[["out_dir"]],
                        paste0("top",params[["top_n"]],"_features_",factor,"_",view,".tsv")),
              sep="\t"
             )
}