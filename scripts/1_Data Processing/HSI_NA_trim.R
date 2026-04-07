
library(terra)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
# Input and output directories
in_dir  <- "E:/Vegetaion Indices Images"
out_dir <- "Spectral_Diversity/VegIndex_NA_trimmed"

#  List all files in directory
allfiles <- list.files(in_dir, full.names = TRUE)

# Exclude sidecar / metadata files
hsi_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]

# Loop through files
for (f in hsi_files) {
  
  message("Processing: ", basename(f))
  
  # Read hyperspectral image
  hsi <- rast(f)
  
  # # Create mask: TRUE where all bands == 0
  # zero_mask <- app(hsi, fun = function(x) {
  #   all(x == 0)
  # })
  
  # Apply mask (convert zero-only pixels to NA) ##((Do this for raw spectra))##
  #hsi_clean <- mask(hsi, zero_mask, maskvalues = 1, updatevalue = NA)
  
  # Count NA values across layers per pixel ##((Do this when dealing with Veg Indices))v####
  na_count <- app(hsi, fun = function(x) sum(is.na(x)))
  # Create mask
  na_mask <- na_count >= 50
  # Apply mask
  hsi_clean <- mask(hsi, na_mask, maskvalues = 1, updatevalue = NA)
  ################################################################## vegIndex Above
  
  
  # Output path (preserve filename)
  out_file <- file.path(out_dir, basename(f))
  
  # Write raster
  writeRaster(
    hsi_clean,
    out_file,
    filetype = "ENVI",
    overwrite = TRUE
  )
  beep()
}

message("All HSI files processed successfully.")
beep(3)

