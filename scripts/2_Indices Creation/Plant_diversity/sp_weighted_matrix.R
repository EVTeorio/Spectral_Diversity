
library(sf)
library(tidyr)
library(dplyr)
library(beepr)

#  SETTINGS
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

crs_proj <- 26916
csv_path      <- "PR_tree_DL.csv"
output_folder <- "Indices_SHPs/"

quads <- st_read("Quad_Scale_SHPs/PR_20m.shp") 
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)


# DATA
tree_df <- read_csv(
  csv_path,
  col_types = cols()              
)

# Keep only rows where cluster_status is "A" or "R"
tree_df <- filter(tree_df, DBH.2024 >= 200 & crown.position %in% c(3, 4, 5))%>%
  filter(cluster_status %in% c("A", "R"))

# POINT LAYER & CROWN BUFFERS
points_sf <- st_as_sf(
  tree_df,
  coords = c("UTMX_CURRENT", "UTMY_CURRENT"), #Changed coord origin
  crs    = crs_proj,
  remove = FALSE
)

points_sf <- points_sf %>%
  mutate(radius_m = cw_m_2025 / 2) %>%
  st_buffer(dist = .$radius_m, endCapStyle = "ROUND")

points_clean <- points_sf %>%
  dplyr::select(sp, quadrat, tag, DBH.2024, crown.position, 
         UTMX_CURRENT, UTMY_CURRENT, radius_m) %>%
  filter(!is.na(UTMX_CURRENT), !is.na(UTMY_CURRENT), !is.na(radius_m)) 
#Go to interactive script for visual##############################

intersections <- st_intersection(points_clean, quads)

intersections <- intersections %>%
  mutate(
    intersect_area = st_area(.)
  )

# Total crown area (original polygons)
crown_areas <- points_clean %>%
  mutate(total_area = st_area(.)) %>%
  st_drop_geometry() %>%
  dplyr::select(tag, total_area)

intersections <- intersections %>%
  left_join(crown_areas, by = "tag") %>%
  mutate(
    prop_area = as.numeric(intersect_area / total_area)
  )

quad_species <- intersections %>%
  st_drop_geometry() %>%
  group_by(Name, sp) %>%
  summarise(
    prop_sum = sum(prop_area, na.rm = TRUE),
    .groups = "drop"
  )

quad_matrix <- quad_species %>%
  pivot_wider(
    names_from = sp,
    values_from = prop_sum,
    values_fill = 0
  )

# Join matrix to quadrat polygons
quads_with_species <- quads %>%
  left_join(quad_matrix, by = "Name")

st_write(
  quads_with_species,
  paste0(output_folder, "species_wieghted_matrix_20m.shp")
)
