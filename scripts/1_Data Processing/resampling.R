
user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

library(terra)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity")

# ---- SETTINGS ----
scale_dirs <- c("10m_smooth", "20m_smooth", "50m_smooth")
scale_cores <- c("10m_smooth" = 8, "20m_smooth" = 4, "50m_smooth" = 2)
quad_spectra_dir <- "Quad_Spectra"
sidecar_pattern <- "\\.(hdr|aux|xml|enp|sta)$"
target_bands <- seq(398, 999, 5)
progress_log <- "logs/resampling_progress.log"

get_band_wavelengths <- function(raster_names) {
  as.numeric(gsub("[^0-9.]", "", raster_names))
}

for (scale_dir in scale_dirs) {
  input_dir <- file.path(quad_spectra_dir, scale_dir)
  output_dir <- file.path(quad_spectra_dir, paste0(scale_dir, "_5nm"))
  n_cores <- min(scale_cores[[scale_dir]], max(1, parallel::detectCores() - 3))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ras_files <- list.files(input_dir, full.names = TRUE)
  tile_files <- ras_files[!grepl(sidecar_pattern, ras_files, ignore.case = TRUE)]
  
  start_message <- paste(
    Sys.time(),
    "Starting resampling for",
    scale_dir,
    "with",
    length(tile_files),
    "rasters using",
    n_cores,
    "cores"
  )
  cat(start_message, "\n")
  cat(start_message, "\n", file = progress_log, append = TRUE)
  
  for (i in seq_along(tile_files)) {
    file <- tile_files[[i]]
    output_file <- file.path(output_dir, basename(file))
    output_hdr <- paste0(output_file, ".hdr")
    
    if (file.exists(output_file) && file.exists(output_hdr)) {
      cat("Skipping existing:", file.path(basename(output_dir), basename(file)), "\n")
      next
    }
    
    tile <- rast(file)
    source_bands <- get_band_wavelengths(names(tile))
    
    if (length(source_bands) != nlyr(tile) || any(is.na(source_bands))) {
      stop("Could not parse source wavelengths for ", file)
    }
    
    resample_pixel <- local({
      source_bands_local <- source_bands
      target_bands_local <- target_bands
      
      function(x) {
        if (all(is.na(x))) {
          return(rep(NA_real_, length(target_bands_local)))
        }
        
        x[x > 1] <- 1.0
        x[x < 0] <- 0.0
        x[is.nan(x) | is.infinite(x)] <- NA_real_
        
        approx(
          x = source_bands_local,
          y = x,
          xout = target_bands_local,
          rule = 1,
          ties = mean
        )$y
      }
    })
    
    resampled_raster <- app(tile, fun = resample_pixel, cores = n_cores)
    names(resampled_raster) <- paste0(target_bands, " nm")
    
    writeRaster(resampled_raster, output_file, filetype = "ENVI", overwrite = TRUE)
    
    progress_message <- paste(
      Sys.time(),
      "Resampled",
      i,
      "of",
      length(tile_files),
      file.path(scale_dir, basename(file))
    )
    cat(progress_message, "\n")
    cat(progress_message, "\n", file = progress_log, append = TRUE)
  }
  
  done_message <- paste(Sys.time(), "Finished resampling for", scale_dir)
  cat(done_message, "\n")
  cat(done_message, "\n", file = progress_log, append = TRUE)
  beep()
}

beep(3)

