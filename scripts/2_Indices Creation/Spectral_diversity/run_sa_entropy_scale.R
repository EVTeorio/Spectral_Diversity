Sys.setenv(RUN_SA_ENTROPY_WORKFLOW = "false")

source("scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R")

selected_scales <- commandArgs(trailingOnly = TRUE)
if (length(selected_scales) == 0) {
  selected_scales <- vapply(SCALE_CONFIG, `[[`, character(1), "scale")
}

selected_config <- Filter(
  function(config) config$scale %in% selected_scales,
  SCALE_CONFIG
)

if (length(selected_config) == 0) {
  stop("No matching scales found. Requested: ", paste(selected_scales, collapse = ", "))
}

run_sa_entropy_workflow(selected_config)
