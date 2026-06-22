Sys.setenv(
  RUN_SA_ENTROPY_WORKFLOW = "false",
  SA_N_BOOT = "2",
  SA_N_SAMPLE = "100",
  SA_N_PAIR_SAMPLE = "1000",
  SA_MAX_EXACT_PAIRS = "10000"
)

source("scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R")

for (config in SCALE_CONFIG) {
  raster_file <- list_raster_files(config$spec_dir)[1]
  cat("\n---", config$scale, basename(raster_file), "---\n")
  result <- process_raster(raster_file, config$scale)
  print(result$summary)
}
