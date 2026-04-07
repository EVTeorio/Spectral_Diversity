
library(terra)
library(sf)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

raster_path <- "Indices_SHPs/Diversity_Rasters/richness_diversity_20m.tif"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_path <- "Indices_SHPs/Diversity_SHPs/richness_diversity_20m.shp"

# READ DATA
quads <- st_read(shp_path)
quads <- quads %>% dplyr::select(-matches("^Dscrptn")) 
div_raster <- rast(raster_path)

# REPROJECT RASTER TO MATCH POLYGON CRS
div_raster_proj <- project(div_raster, crs(vect(quads)))

# EXTRACT MEAN PIXEL VALUE PER POLYGON
quads_vect <- vect(quads)
extracted  <- terra::extract(div_raster_proj, quads_vect, fun = mean, na.rm = TRUE)

# ADD TO SHAPEFILE AS NEW COLUMN
quads <- quads %>%
  mutate(richness_diversity = extracted[, 2])

# WRITE OUT
st_write(quads, out_path, delete_layer = TRUE)

