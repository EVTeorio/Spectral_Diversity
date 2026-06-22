user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

library(terra)
library(signal)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

# --- Parameters ---
sg_p <- 3   # polynomial order
sg_n <- 15  # window size

scale_dirs <- c("10m", "20m", "50m")
scale_cores <- c("10m" = 8, "20m" = 4, "50m" = 2)
quad_spectra_dir <- "Spectral_Diversity/Quad_Spectra"
sidecar_pattern <- "\\.(hdr|aux|xml|enp|sta)$"
progress_log <- "Spectral_Diversity/logs/smoothing_progress.log"

smooth_spectrum <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  
  signal::sgolayfilt(x, p = 3, n = 15)
}

for (scale_dir in scale_dirs) {
  in_dir <- file.path(quad_spectra_dir, scale_dir)
  out_dir <- file.path(quad_spectra_dir, paste0(scale_dir, "_smooth"))
  n_cores <- min(scale_cores[[scale_dir]], max(1, parallel::detectCores() - 3))
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  ras_files <- list.files(in_dir, full.names = TRUE)
  raster_files <- ras_files[!grepl(sidecar_pattern, ras_files, ignore.case = TRUE)]
  
  start_message <- paste(
    Sys.time(),
    "Starting smoothing for",
    scale_dir,
    "with",
    length(raster_files),
    "rasters using",
    n_cores,
    "cores"
  )
  cat(start_message, "\n")
  cat(start_message, "\n", file = progress_log, append = TRUE)
  
  for (i in seq_along(raster_files)) {
    f <- raster_files[[i]]
    out_file <- file.path(out_dir, basename(f))
    out_hdr <- paste0(out_file, ".hdr")
    
    if (file.exists(out_file) && file.exists(out_hdr)) {
      cat("Skipping existing:", file.path(basename(out_dir), basename(f)), "\n")
      next
    }
    
    r <- rast(f)
    orig_names <- names(r)
    
    r_smoothed <- app(r, fun = smooth_spectrum, cores = n_cores)
    names(r_smoothed) <- orig_names
    
    writeRaster(r_smoothed, out_file, filetype = "ENVI", overwrite = TRUE)
    progress_message <- paste(Sys.time(), "Smoothed", i, "of", length(raster_files), file.path(scale_dir, basename(f)))
    cat(progress_message, "\n")
    cat(progress_message, "\n", file = progress_log, append = TRUE)
  }
  
  done_message <- paste(Sys.time(), "Finished smoothing for", scale_dir)
  cat(done_message, "\n")
  cat(done_message, "\n", file = progress_log, append = TRUE)
  beep()
}

beep(3)
