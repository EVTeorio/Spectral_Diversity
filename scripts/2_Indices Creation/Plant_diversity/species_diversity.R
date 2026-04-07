
library(sf)
library(dplyr)

# --- SETTINGS ---
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

shapefile_path <- "Indices_SHPs/species_wieghted_matrix_20m.shp"
crs_proj <- 26916
# --- READ DATA ---
quads <- st_read("Quad_Scale_SHPs/PR_20m.shp")
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)
quads_sf <- st_read(shapefile_path)

# List of species columns (adjust based on your data)
species_cols <- quads_sf %>%
  dplyr::select(-Name) %>%
  st_drop_geometry() %>%
  colnames()

# --- FUNCTION TO CALCULATE DIVERSITY METRICS ---
calculate_diversity <- function(abundances) {
  # Remove zero abundances
  abundances <- abundances[abundances > 0]
  
  if(length(abundances) == 0) return(
    data.frame(richness = 0, shannon = 0, simpson = 0, evenness = 0)
  )
  
  # Relative abundances
  p <- abundances / sum(abundances)
  
  richness <- length(abundances)
  shannon  <- -sum(p * log(p))
  simpson  <- 1 - sum(p^2)
  evenness <- ifelse(richness > 1, shannon / log(richness), 0)
  
  data.frame(richness, shannon, simpson, evenness)
}

# --- APPLY TO EACH POLYGON ---
diversity_df <- quads_sf %>%
  rowwise() %>%
  mutate(
    metrics = list(calculate_diversity(c_across(all_of(species_cols))))
  ) %>%
  tidyr::unnest(metrics)%>%
  st_drop_geometry()

# Optional: attach metrics back to shapefile
quads_final <- quads %>%
  bind_cols(diversity_df %>% dplyr::select(richness, shannon, simpson, evenness))

st_write(quads_final, "Indices_SHPs/Diversity_SHPs/species_diversity_20m.shp")



