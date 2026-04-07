
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(raster)
library(tidyr)
library(terra)

#  SETTINGS
library(beepr)
csv_path      <- "C:/Users/PaintRock/Downloads/PR_tree_DL.csv"
output_folder <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/Test"

crs_proj <- 26916

diversity_index <- "simpson"
window_radius_m <- 20


# #  SETTINGS
# csv_path      <- "C:/PRFDP/recensus/data/PR_tree_DL.csv"
# output_folder <- "C:/PRFDP/recensus/output"
# dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
# 
# 
# crs_proj <- 26916 # changed CRS
# 
# diversity_index <- "simpson"
# window_radius_m <- 20
# setwd("C:/PRFDP/recensus/data/analysis")

# DATA

tree_df <- read_csv(
  csv_path,
  col_types = cols()              
)


# Keep only rows where cluster_status is "A" or "R"
tree_df <- tree_df %>%
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

# BUILD 1 × 1 m

bbox <- st_bbox(points_sf)

grid_1m <- st_make_grid(
  x        = st_as_sfc(bbox),
  cellsize = 1,         
  square   = TRUE
)

grid_sf <- st_sf(
  geometry = grid_1m,
  crs      = st_crs(points_sf)
) %>%
  mutate(cell_id = row_number())  

# moving window

centroids <- st_centroid(grid_sf)
window_sf <- st_buffer(centroids, dist = window_radius_m)

idx <- st_intersects(window_sf, points_sf)


# DIVERSITY METRIC
species_vec <- points_sf$sp

diversity_vec <- vapply(
  idx,
  function(i) {
    if (length(i) == 0) return(0)
    
    sp_counts <- table(species_vec[i])
    p <- sp_counts / sum(sp_counts)
    
    switch(
      diversity_index,
      
      "richness" = length(sp_counts),
      
      "shannon"  = -sum(p * log(p)),
      
      "simpson"  = 1 - sum(p^2),
      
      "evenness" = {
        H <- -sum(p * log(p))
        H / log(length(sp_counts))
      },
      
      stop("Unknown diversity_index: ", diversity_index)
    )
  },
  numeric(1)
)

#Grid
col_name <- paste0(diversity_index, "_", window_radius_m, "m")
grid_sf[[col_name]] <- diversity_vec

grid_vect <- vect(grid_sf)

r_template <- rast(
  ext(grid_vect),
  resolution = 1,
  crs = crs(grid_vect)
)


diversity_rast <- rasterize(
  grid_vect,
  r_template,
  field = col_name
)

out_raster <- paste0(
  diversity_index,
  "_diversity_",
  window_radius_m,
  "m.tif"
)

writeRaster(
  diversity_rast,
  filename = out_raster,
  overwrite = TRUE
)

# Map


if (diversity_index == "shannon") {
  fill_limits <- c(0, 3.5)
  fill_colors <- viridis::viridis(5, option = "C")
} else if (diversity_index %in% c("simpson", "evenness")) {
  fill_limits <- c(0, 1)
  fill_colors <- viridis::viridis(5, option = "C")
} else {
  fill_limits <- range(grid_sf[[col_name]], na.rm = TRUE)
  fill_colors <- viridis::viridis(5)
}


diversity_map <- ggplot() +
  geom_sf(
    data = grid_sf,
    aes(fill = .data[[col_name]]),
    color = NA
  ) +
  scale_fill_gradientn(
    colours = fill_colors,
    limits  = fill_limits,
    oob     = scales::squish,
    name    = paste0(
      tools::toTitleCase(diversity_index),
      " diversity\n(",
      window_radius_m,
      " m neighbourhood)"
    ),
    guide   = guide_colourbar(
      barwidth  = 0.5,
      barheight = 15
    )
  ) +
  coord_sf(crs = st_crs(grid_sf)) +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    axis.title      = element_blank()
  ) +
  ggtitle(
    paste0(
      tools::toTitleCase(diversity_index),
      " Diversity (",
      window_radius_m,
      " m Moving Window, 1 m Grid)"
    )
  )

print(diversity_map)

