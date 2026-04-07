
setwd("E:/Git Paint Rock 1.0/lecospec")
source("Functions/lecospectR.R")

library(terra)
library(spectrolab)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

  # ---- SETTINGS ----
input_dir  <- "Spectral_Diversity/Quad_Spectra/10m"
output_dir <- "Spectral_Diversity/Quad_Spectra/10m_resampled_5nm"

# List raster tiles
ras_files <- list.files(input_dir, full.names = TRUE)
tile_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

# ---- PROCESS LOOP ----
for (file in tile_files) {
  
  cat("Processing:", basename(file), "\n")
  
  tile <- rast(file)
  
  # Read tile
  #tile <- writeRaster(tile, tempfile(fileext = ".tif"), overwrite = TRUE)
  
  # Convert raster to dataframe (pixels x bands)
  df <- as.data.frame(tile, na.rm = FALSE)
  df <- filter_bands(df)
  df <- df_to_speclib(df, type = "spectrolab")
  
  # Resample
  df_resampled <- spectrolab::resample(
    df,
    new_bands = seq(398, 999, 5),
    fwhm = 1
  )
  
  # Convert back to matrix
  mat_resampled <- as.matrix(df_resampled)
  
  # Create empty raster with same geometry
  resampled_raster <- rast(
    nrows = nrow(tile),
    ncols = ncol(tile),
    nlyrs = ncol(mat_resampled),
    crs   = crs(tile),
    ext   = ext(tile)
  )

  values(resampled_raster) <- mat_resampled

  names(resampled_raster) <- paste0(seq(398, 999, 5), " nm")

  output_file <- file.path(output_dir, basename(file))
  
  # Write to disk
  writeRaster(resampled_raster, output_file, filetype="ENVI", overwrite = TRUE)
}
beep()

