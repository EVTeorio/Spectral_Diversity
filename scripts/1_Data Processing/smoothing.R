library(terra)
library(signal)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

# --- Paths ---
in_dir  <- "Spectral_Diversity/Quad_Spectra/20m_resampled_5nm"
out_dir <- "Spectral_Diversity/Quad_Spectra/20m_smoothed_5nm"

# --- Parameters ---
sg_p <- 3   # polynomial order
sg_n <- 7  # window size

# --- List files ---
ras_files <- list.files(in_dir, full.names = TRUE)
raster_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

# --- Loop through files ---
for(f in raster_files){
  r <- rast(f)
  
  # Store original band names
  orig_names <- names(r)
  
  # Apply Savitzky-Golay smoothing band by band
  r_smoothed <- app(r, fun = function(x){
    if(all(is.na(x))){
      return(x)  # keep NA if all values are NA
    } else {
      return(sgolayfilt(x, p = sg_p, n = sg_n))
    }
  })
  
  # Restore original band names
  names(r_smoothed) <- orig_names
  
  # --- Save output ---
  out_file <- file.path(out_dir, basename(f))
  writeRaster(r_smoothed, out_file, filetype = "ENVI")
  cat("Smoothed:", basename(f), "\n")
  beep()
}

beep(3)
