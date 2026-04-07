

library(terra)

HSI <- rast("E:/HSI_Files_Parsing/raw_11159_rd_rf_or")

# HSI is your SpatRaster
r <- HSI

# Get current extent
e <- ext(r)

# Total width in x-direction
x_width <- e$xmax - e$xmin

# Central 40% bounds
x_min_new <- e$xmin + 0.30 * x_width
x_max_new <- e$xmin + 0.70 * x_width


# If pixel based is needed   ########

# ncols <- ncol(r)
# 
# start_col <- floor(ncols * 0.30)
# end_col   <- ceiling(ncols * 0.70)
# 
# HSI_mid40 <- r[, start_col:end_col]
###########################################

# New extent (keep full y-range)
e_mid40 <- ext(
  x_min_new,
  x_max_new,
  e$ymin,
  e$ymax
)

# Crop raster
HSI_mid40 <- crop(r, e_mid40)

writeRaster(HSI_mid40, "SOFOR/HSI_40%_cropped/raw_11159_rd_rf_or.envi")


########################################################################3
library(terra)

r <- HSI

# Original extent
e <- ext(r)

# Center of the raster
x_center <- (e$xmin + e$xmax) / 2
y_center <- (e$ymin + e$ymax) / 2

# Convert meters to degrees longitude at this latitude
meters_per_degree_lon <- 111320 * cos(y_center * pi / 180)

# Half-width in degrees (50 m total → 25 m each side)
half_width_deg <- 30 / meters_per_degree_lon

# New extent: middle 50 meters, full height
e_50m <- ext(
  x_center - half_width_deg,
  x_center + half_width_deg,
  e$ymin,
  e$ymax
)

# Crop
HSI_mid50m <- crop(r, e_50m)

plot(HSI_mid50m[[328]])

writeRaster(HSI_mid50m, "SOFOR/HSI_40%_cropped/raw_11159_rd_rf_or.envi", )

