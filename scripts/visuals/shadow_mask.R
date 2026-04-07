
library(terra)
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")


hsi <- rast("Quad_Spectra/20m_smoothed_5nm/604")
rgb <- rast("Quad_Spectra/20m_RGB/604.tif")

# Parameters
# -----------------------
best_band   <- "563 nm"
threshold   <- 0.0305476
direction   <- "<"

# -----------------------
# Select Band
# -----------------------
band_603 <- hsi[[best_band]]

# -----------------------
# Apply Classification Rule
# -----------------------
if (direction == ">") {
  
  # Lower values = Sunlit
  classified <- ifel(band_603 < threshold, 1, 2)
  
} else {
  
  # Higher values = Sunlit
  classified <- ifel(band_603 > threshold, 1, 2)
}

plot(classified)

# Assign class labels
levels(classified) <- data.frame(
  value = c(1, 2),
  class = c("Sunlit", "Shadow")
)

# Rename layer
names(classified) <- "Illumination_Class"

# -----------------------
# Optional: Create Mask (Sunlit Only)
# -----------------------
sunlit_mask <- classified == 1
plot(sunlit_mask)

# Apply mask to RGB
rgb_sunlit <- mask(rgb, sunlit_mask, maskvalues = 0)
plot(rgb_sunlit)

plotRGB(rgb, r = 1, g = 2, b = 3, stretch = "lin")
plotRGB(rgb_sunlit, r = 1, g = 2, b = 3, stretch = "lin")

# Set plotting layout to 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))

# Original RGB
plotRGB(rgb,
        r = 1, g = 2, b = 3,
        stretch = "lin",
        main = "Original RGB")

# Sunlit Masked RGB
plotRGB(rgb_sunlit,
        r = 1, g = 2, b = 3,
        stretch = "lin",
        main = "Sunlit Only")

# Reset layout
par(mfrow = c(1, 1))
# -----------------------
# Write Outputs
# -----------------------

# Full classification raster
writeRaster(classified,
            "hsi_illumination_class.tif",
            overwrite = TRUE)
