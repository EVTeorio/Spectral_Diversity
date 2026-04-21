
library(raster)
library(sf)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")
crs_proj <- 26916
# 1) Read in the ADF file as a raster
raster_file <- "Indices_SHPs/Other_variables/lci_update/w001001x.adf"
raster_layer <- raster(raster_file)

# 2) Read in the shapefile (quads) as a Simple Feature (SF) object
quads <- st_read("Quad_Scale_SHPs/PR_20m.shp")  # Update with your shapefile path
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)

# Extract raster values within each quad polygon
extracted_values <- extract(raster_layer, quads, df = FALSE)  

# 2) For each quad, calculate the dominant (most frequent) categorical value
dominant_values <- lapply(extracted_values, function(vals) {
  if (length(vals) > 0) {
    # Create a frequency table of the pixel values
    value_counts <- table(vals)  # Each unique pixel value is a row, and frequency is counted
    
    # Find the pixel value with the highest frequency (most common)
    dominant_value <- names(value_counts)[which.max(value_counts)]
    
    return(dominant_value)
  } else {
    # If no values are found (e.g., quad doesn't overlap with any raster pixels), return NA
    return(NA)
  }
})

# 3) Add the dominant values as a new column in the 'quads' SF object
quads$dominant_value <- unlist(dominant_values)

# 4) View the result
head(quads)

# 5) Optionally, save the result to a new shapefile
st_write(quads, "Indices_SHPs/Other_variables/20mTopo_Index.shp")
