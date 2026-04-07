
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(vegan)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
justforcrs <- rast("E:/Updated LiDAR/PRFPD_CHM_leafOff.tiff")
#--------------------------------------------------
# 1. Read and clean quad polygons
#--------------------------------------------------
quads <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp")

quads_clean <- quads %>% select(-matches("^Dscrptn")) #%>%
 # mutate(Name = sub_id) %>%   # replace Name with sub_id
 # select(-sub_id)

#--------------------------------------------------
# 2. Read census table and convert to sf points
#--------------------------------------------------
census <- read.csv("Spectral_Diversity/PR_tree_DL.csv")

# Filter which tree Size and or CPI
census <- filter(census, DBH.2024 >= 200 | crown.position %in% c(4, 5))

census_sf <- st_as_sf(
  census,
  coords = c("UTMX_CURRENT", "UTMY_CURRENT"),
  crs = st_crs(justforcrs)   # NAD83 / UTM zone 16N
)
census_sf <- st_transform(census_sf, st_crs(quads))
crs(census_sf)

#--------------------------------------------------
# 3. Spatially join points to quads
#    (each tree gets the quad it falls in)
#--------------------------------------------------
census_joined <- st_join(
  census_sf,
  quads_clean,
  join = st_within
)

table(is.na(census_joined$Join_ID))

#--------------------------------------------------
# 4. Species counts per quad (spatially derived)
#--------------------------------------------------
species_counts <- census_joined %>%
  st_drop_geometry() %>%
  filter(!is.na(Name)) %>%
  group_by(Name, sp) %>%
  summarise(n = n(), .groups = "drop")

species_wide <- species_counts %>%
  pivot_wider(
    names_from = sp,
    values_from = n,
    values_fill = 0
  )

#--------------------------------------------------
# 5. Join species matrix back to quads
#--------------------------------------------------
quads_species <- quads_clean %>%
  left_join(species_wide, by = "Name")

#--------------------------------------------------
# 6. Diversity indices
#--------------------------------------------------
species_matrix <- quads_species %>%
  st_drop_geometry() %>%
  select(where(is.numeric), -Name) %>%
  as.matrix()

quads_species$shannon <- diversity(species_matrix, index = "shannon")
quads_species$simpson <- diversity(species_matrix, index = "simpson")

#--------------------------------------------------
# 7. Write output
#--------------------------------------------------
st_write(quads_species, "Spectral_Diversity/Indices_SHPs/sp_diversity_20m_200DBH.shp", delete_layer = TRUE)

##############################################################











