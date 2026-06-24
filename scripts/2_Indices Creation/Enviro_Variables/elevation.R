

library(raster)
library(sf)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")
crs_proj <- 26916

# Set the file paths for the DTM raster and the shapefile (quads)
dtm_file <- "E:/Updated LiDAR/PRFPD_DTM_leafOff.tiff"

#Read in the DTM raster file (Digital Terrain Model)
dtm_raster <- raster(dtm_file)

#Read in the shapefile (quads)
quads <- st_read("Quad_Scale_SHPs/PR_20m.shp")  # Update with your shapefile path
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)


# Extracting elevation values for each quad (polygon)
extracted_values <- extract(dtm_raster, quads, df = FALSE) 

# For each quad, calculate the average elevation and the elevation range (max - min)
elevation_stats <- lapply(extracted_values, function(vals) {
  if (length(vals) > 0) {
    # Calculate the average elevation
    avg_elevation <- mean(vals, na.rm = TRUE)
    
    # Calculate the elevation range (max - min)
    elevation_range <- max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)
    
    return(c(avg_elevation = avg_elevation, elevation_range = elevation_range))
  } else {
    # If no values are found (quad doesn't overlap with valid pixels), return NA
    return(c(avg_elevation = NA, elevation_range = NA))
  }
})

# Convert the list of elevation stats into a data frame and add it to the quads shapefile
elevation_stats_df <- do.call(rbind, elevation_stats) 
elevation_stats_df <- as.data.frame(elevation_stats_df)
colnames(elevation_stats_df) <- c("avg_elevation", "elevation_range")

# Add the calculated stats as new columns in the quads SF object
quads <- cbind(quads, elevation_stats_df)
head(quads)

# Optionally, save the updated shapefile with the new statistics
st_write(quads, "Indices_SHPs/Other_variables/20m_elevation.shp")
