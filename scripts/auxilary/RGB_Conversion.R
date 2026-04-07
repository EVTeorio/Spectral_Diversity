

setwd("E:/Git Paint Rock 1.0/lecospec")
source("Functions/lecospectR.R")

library(terra)
library(spectrolab)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

# ---- SETTINGS ----
input_dir  <- "Spectral_Diversity/Quad_Spectra/50m_resampled_5nm"
output_dir <- "Spectral_Diversity/Quad_Spectra/50m_RGB"

# List raster tiles
ras_files <- list.files(input_dir, full.names = TRUE)
tile_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

# Loop over each raster tile
lapply(tile_files, function(file) {
  
  cat("Processing file:", basename(file), "\n")
  
  # Load raster
  HSI <- rast(file)
  
  # Select bands
  red_band   <- 130
  green_band <- 84
  blue_band  <- 22
  
  RGB <- HSI[[c(red_band, green_band, blue_band)]]
  
  normalize_layer <- function(x) {
    rmin <- global(x, "min", na.rm=TRUE)[1,1]
    rmax <- global(x, "max", na.rm=TRUE)[1,1]
    (x - rmin) / (rmax - rmin)
  }
  
  RGB_norm <- c(
    normalize_layer(RGB[[1]]),
    normalize_layer(RGB[[2]]),
    normalize_layer(RGB[[3]])
  )
  names(RGB_norm) <- c("Red", "Green", "Blue")
  
  # Construct output filename
  output_name <- paste0(tools::file_path_sans_ext(basename(file)), ".tif")
  output_path <- file.path(output_dir, output_name)
  
  # Write the RGB raster
  writeRaster(RGB_norm, output_path, overwrite = TRUE)
  
})

beep()

HS <- rast(tile_files[[6]])
plot(HS)
ncell(HS)
global(HS[[1]], fun="notNA")
