library(UpSetR)
library(jsonlite)
library(yaml)
#library(ggplot2)
#library(tidyverse)

build_combinedModuleList<-function(stages){
    #Parameters
    # stages:named list of modules containing filepath to model module json file containing gene/protein names/symbols
    combined_modules<-c()
    for (stage in names(stages)){
        modules_in_stage<-read_json(stages[[stage]], simplifyVector = TRUE)
        for (module in names(modules_in_stage)){
            new_module<-paste0(stage,'_',module)
            combined_modules[new_module]<-modules_in_stage[module]
        }
    }
    combined_modules
}

### MAIN CODE ###

#read dvc parameters YAML
params = yaml::read_yaml("params.yaml")[["WGCNA"]][["upset"]]

stages=params[["stages"]]
print(stages)

#built list containing all modules of all mouse models
modules<-fromList(build_combinedModuleList(stages))
#build upset plots
dir.create(params[["out_dir"]])
pdf(file.path(params[["out_dir"]],"ModuleUpset.pdf"))
upset(modules,
      nsets=dim(modules)[2],
      order.by = params[["order_by"]],
      nintersects= params[["nintersects"]], 
      mb.ratio = c(params[["matrix_plot_ratio"]], 
                   params[["main_bar_plot_ratio"]]),
      text.scale = params[["text_scale"]])
dev.off()
png(file.path(params[["out_dir"]],"ModuleUpset.png"),
    width = params[["png_width"]], 
    height = params[["png_height"]],
    res=params[["png_res"]])
upset(modules,nsets=dim(modules)[2],
      order.by = params[["order_by"]],
      nintersects=params[["nintersects"]], 
      mb.ratio = c(params[["matrix_plot_ratio"]], 
                   params[["main_bar_plot_ratio"]]),
      text.scale = params[["text_scale"]])
dev.off()