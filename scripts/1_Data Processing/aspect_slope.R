
library(terra)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

dtm <- rast("E:/Updated LiDAR/PRFPD_DTM_leafOff.tiff")
dtm

# Compute slope in radians first
slope_rad <- terrain(dtm, v = "slope", unit = "radians", neighbors = 8)
beep()
# Convert to percent rise
slope_percent <- tan(slope_rad) * 100
plot(slope_rad)

aspect_deg <- terrain(dtm, v = "aspect", unit = "degrees", neighbors = 8)
aspect_rad <- aspect_deg * pi / 180
aspect_sin <- sin(aspect_rad)  # east-west structure
aspect_cos <- cos(aspect_rad)  # north-south structure
plot(aspect_sin)

#Topographic Exposure Index
slope_tan <- tan(slope_rad)
tei_north <- aspect_cos * slope_tan
tei_east <- aspect_sin * slope_tan
plot(tei_east)

# exposure_mag <- sqrt(tei_east^2 + tei_north^2)
# plot(exposure_mag)


writeRaster(slope_percent,
            "Quad_Values/Other_variables/Slope_Aspect_rasters/slope_percent.tif",
            overwrite = TRUE)

writeRaster(aspect_cos,
            "Quad_Values/Other_variables/Slope_Aspect_rasters/cos_aspect.tif",
            overwrite = TRUE)

writeRaster(aspect_sin,
            "Quad_Values/Other_variables/Slope_Aspect_rasters/sin_aspect.tif",
            overwrite = TRUE)

writeRaster(tei_north,
            "Quad_Values/Other_variables/Slope_Aspect_rasters/tei_north.tif",
            overwrite = TRUE)

writeRaster(tei_east,
            "Quad_Values/Other_variables/Slope_Aspect_rasters/tei_east.tif",
            overwrite = TRUE)

###############################################################################
#######   Breaking into Quads #################################################
###############################################################################

read.csv()

crs_proj <- 26916
quads <- st_read("Quad_Scale_SHPs/PR_20m.shp")
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)

topo_stack <- c(
  slope_percent,
  aspect_sin,
  aspect_cos,
  tei_north,
  tei_east
)

names(topo_stack) <- c(
  "slope_percent",
  "aspect_sin",
  "aspect_cos",
  "tei_north",
  "tei_east"
)

# -----------------------------
# Zonal extraction (mean per quad)
# -----------------------------
z <- terra::extract(
  topo_stack,
  quads,
  fun = mean,
  na.rm = TRUE,
  cells = FALSE
)

df <- data.frame(
  name = quads$Name,
  z[, -1]   # remove ID column from extract output
)

write.csv(df,"Quad_Values/Other_variables/Slope_Aspect.csv" )
