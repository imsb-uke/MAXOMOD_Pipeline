source(here::here("src/shared/R-helpers/get_args.R"))
source(here::here("src/shared/R-helpers/update_list.R"))

suppressPackageStartupMessages({
  requireNamespace("yaml")
})


get_settings <- function(omic, stage, dataset_name=NULL) {
  if (is.null(dataset_name)) {
    cmdline_args <- get_args()
    dataset_name <- cmdline_args[["dataset_name"]]
  }

  params <- yaml::read_yaml("params.yaml")
  params_for_stage <- params[[omic]][[stage]]

  our_settings <- list()

  def_settings <- params_for_stage[["default_settings"]]
  our_settings <- update_list(our_settings, def_settings)

  dataset_settings <- params_for_stage[["stages"]][[dataset_name]]
  if ("settings" %in% names(dataset_settings)) {
    our_settings <- update_list(our_settings, dataset_settings[["settings"]])
  }

  def_syssettings <- params_for_stage[["default_system_settings"]]
  our_settings <- update_list(our_settings, def_syssettings)

  if ("system_settings" %in% names(dataset_settings)) {
    our_settings <- update_list(our_settings, dataset_settings[["system_settings"]])
  }
  our_settings
}
