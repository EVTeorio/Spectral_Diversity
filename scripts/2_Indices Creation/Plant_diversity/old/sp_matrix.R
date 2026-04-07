

library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(terra)

# SETTINGS
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
csv_path  <- "Spectral_Diversity/PR_tree_DL.csv"
shp_path  <- "Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp"
out_path  <- "Spectral_Diversity/Indices_SHPs/Diversity_SHPs/PR_20m_sp_matrix.shp"
crs_proj  <- 26916

# DATA
tree_df <- filter(tree_df, DBH.2024 >= 200 | crown.position %in% c(4, 5))%>%
  filter(cluster_status %in% c("A", "R"))

# CROWN BUFFERS
points_sf <- st_as_sf(
  tree_df,
  coords = c("UTMX_CURRENT", "UTMY_CURRENT"),
  crs    = crs_proj,
  remove = FALSE
) %>%
  mutate(radius_m = cw_m_2025 / 2) %>%
  st_buffer(dist = .$radius_m, endCapStyle = "ROUND")

# READ & CLEAN QUADS
quads <- st_read(shp_path, quiet = TRUE) %>%
  select(-matches("^Dscrptn")) %>%
  st_transform(crs = crs_proj)

# INTERSECT CROWNS WITH QUADS
intersected <- st_join(quads, points_sf, join = st_intersects)

# COUNT SPECIES PER QUAD
species_counts <- intersected %>%
  st_drop_geometry() %>%
  filter(!is.na(sp)) %>%
  group_by(Name, sp) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = sp,
    values_from = n,
    values_fill = 0
  )

# JOIN COUNTS BACK TO QUADS
quads_out <- quads %>%
  left_join(species_counts, by = "Name")

# WRITE OUT
st_write(quads_out, out_path, delete_layer = TRUE)
