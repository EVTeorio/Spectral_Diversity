
library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

# --- Load Rasters ---
resampled <- rast("Quad_Spectra/20m_resampled_5nm/0")
smoothed <- rast("Quad_Spectra/20m_smoothed_5nm/0")

# --- Function to plot spectral signatures for two rasters side by side ---
plot_compare_spectra <- function(r1, r2, n = 1, seed = 123) {
  
  set.seed(seed)  # reproducible pixel selection
  
  # Pick n random valid pixels based on first raster
  valid_cells <- which(!is.na(values(r1[[1]])))
  if(length(valid_cells) < n) stop("Not enough valid pixels!")
  
  pixel_cells <- sample(valid_cells, n)
  
  # Extract values from both rasters
  vals1 <- terra::extract(r1, pixel_cells)[,-1]
  vals2 <- terra::extract(r2, pixel_cells)[,-1]
  
  # Extract wavelengths from column names (use r1)
  wavelengths <- as.numeric(gsub("[^0-9.]", "", colnames(vals1)))
  
  # Convert to long format for ggplot
  df1 <- as.data.frame(vals1)
  df1$Pixel <- paste0("Pixel_", 1:nrow(df1))
  df1$Raster <- "Smoothed"
  
  df2 <- as.data.frame(vals2)
  df2$Pixel <- paste0("Pixel_", 1:nrow(df2))
  df2$Raster <- "Resampled"
  
  df_long <- bind_rows(df1, df2) %>%
    pivot_longer(cols = -c(Pixel, Raster),
                 names_to = "Wavelength",
                 values_to = "Reflectance") %>%
    mutate(Wavelength = as.numeric(gsub("[^0-9.]", "", Wavelength)))
  
  # --- Plot locator map for the first raster only ---
  plot(r1[[1]], main = "Selected Pixel(s)")
  points(xFromCell(r1, pixel_cells), yFromCell(r1, pixel_cells),
         col = "red", pch = 16, cex = 1.5)
  
  # --- ggplot of spectral signatures side by side ---
  ggplot(df_long, aes(x = Wavelength, y = Reflectance, color = Pixel)) +
    geom_line(size = 1) +
    facet_wrap(~Raster) +
    labs(
      title = paste("Comparison of Spectral Signatures"),
      x = "Wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom")
}

# --- Run the comparison ---
plot_compare_spectra(smoothed, resampled, n = 1, seed = 16)

