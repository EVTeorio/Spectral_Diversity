
library(terra)
library(sf)
library(dplyr)
library(snow)
library(beepr)
library(parallel)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

hsi_dir  <- "Spectral_Diversity/Quad_Spectra/20m"
out_dir  <- "Spectral_Diversity/Quad_Spectra/Mosiac/"

# List hyperspectral tiles
# -----------------------------
allfiles <- list.files(hsi_dir, full.names = TRUE)
hsi_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]

# -----------------------------
# Build virtual mosaic (NO resampling)
# -----------------------------
v <- vrt(hsi_files)
beep(3)
# -----------------------------
# Write to ENVI (single pass, on disk)
# -----------------------------
writeRaster(
  v,
  file.path(out_dir, "mosaic_20m"),
  filetype  = "ENVI",
  overwrite = TRUE)

beep(3) #lets me know when it done

# Stop cluster
stopCluster(cl)
