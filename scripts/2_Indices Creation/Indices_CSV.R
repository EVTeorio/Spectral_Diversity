
library(sf)
library(dplyr)
library(tools)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

spectral_dir <- "Spectral_Diversity/Indices_SHPs/Spectral_diversitySHPs"
species_dir  <- "Spectral_Diversity/Indices_SHPs/Diversity_SHPs"
sp_matrix <- st_read("Spectral_Diversity/Indices_SHPs/PR_20m_sp_matrix.shp")
# Drop geometry from sp_matrix
sp_matrix_data <- st_drop_geometry(sp_matrix)

# Get shapefile paths
spectral_files <- list.files(spectral_dir, pattern = "\\.shp$", full.names = TRUE)
species_files  <- list.files(species_dir, pattern = "\\.shp$", full.names = TRUE)

# Read shapefiles
spectral <- lapply(spectral_files, st_read, quiet = TRUE)
species  <- lapply(species_files, st_read, quiet = TRUE)

for (i in seq_along(spectral)) {
  new_name <- file_path_sans_ext(basename(spectral_files[i]))
  names(spectral[[i]])[names(spectral[[i]]) == "spctrl_"] <- new_name
}
spectral <- spectral[-2]



# Function to extract Name + second column as a data.frame
select_data <- function(shp) {
  col_name <- names(shp)[2]
  # Keep Name + second column
  df <- shp %>% st_drop_geometry() %>% dplyr::select(Name, all_of(col_name))
  return(df)
}

# Extract the second columns
spectral_data <- lapply(spectral, select_data)
species_data  <- lapply(species, select_data)

# Combine all into one list
all_data <- c(spectral_data, species_data)

# Merge all data.frames by Name iteratively
merged_data <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), all_data)

# Merge with your previously combined spectral + species dataframe
merged_data <- merged_data %>%
  left_join(sp_matrix_data, by = "Name")

# Attach geometry from the first shapefile
merged_sf <- st_sf(merged_data, geometry = spectral[[1]]$geometry)
# Convert geometry to WKT string
merged_sf$geometry <- st_as_text(merged_sf$geometry, digits = 16)

write.csv(merged_sf, "Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv")
