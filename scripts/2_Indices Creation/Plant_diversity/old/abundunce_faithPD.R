
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(tidyr)
library(terra)
library(beepr)
library(ape)
library(picante)
library(V.PhyloMaker2)


setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

# SETTINGS
csv_path         <- "Spectral_Diversity/PR_tree_DL.csv"
taxa_path        <- "Spectral_Diversity/51sp_taxanomy.csv"
output_folder    <- "Spectral_Diversity/Indices_SHPs/Diversity_Rasters"
crs_proj         <- 26916
window_radius_m  <- 20

# DATA
tree_df <- read_csv(csv_path, col_types = cols())
tree_df <- filter(tree_df, DBH.2024 >= 200 | crown.position %in% c(4, 5))%>%
  filter(cluster_status %in% c("A", "R"))

taxa_table <- read.csv(taxa_path)

# BUILD TIME-CALIBRATED PHYLOGENETIC TREE VIA V.PhyloMaker2
phylo_input <- taxa_table %>%
  mutate(species = paste(genus, species, sep = " ")) %>%
  select(species, genus, family) %>%
  filter(species != "Carya sp") %>%       # Drop unknown sp
  distinct(species, .keep_all = TRUE)     # Drop duplicate binomials

phylo_tree <- phylo.maker(
  sp.list   = phylo_input,
  tree      = GBOTB.extended.TPL,
  nodes     = nodes.info.1.TPL,
  scenarios = "S3"
)$scenario.3
beep()
# Map sp_code to tip labels (genus_species format matching tree tips)
sp_to_tip <- setNames(
  paste(taxa_table$genus, taxa_table$species, sep = "_"),
  taxa_table$sp_code
)
all_tips <- phylo_tree$tip.label

# POINT LAYER & CROWN BUFFERS
points_sf <- st_as_sf(
  tree_df,
  coords = c("UTMX_CURRENT", "UTMY_CURRENT"),
  crs    = crs_proj,
  remove = FALSE
)
points_sf <- points_sf %>%
  mutate(
    radius_m  = cw_m_2025 / 2,
    tip_label = sp_to_tip[sp]
  ) %>%
  filter(!is.na(tip_label)) %>%
  filter(tip_label %in% all_tips) %>%
  st_buffer(dist = .$radius_m, endCapStyle = "ROUND")

# BUILD 1 x 1 m GRID
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

# MOVING WINDOW INTERSECTIONS
centroids <- st_centroid(grid_sf)
window_sf <- st_buffer(centroids, dist = window_radius_m)
idx       <- st_intersects(window_sf, points_sf)
beep()
tip_vec <- points_sf$tip_label

pd_vec <- vapply(
  idx,
  function(i) {
    
    if (length(i) == 0) return(0)
    
    tips <- tip_vec[i]
    tips <- tips[tips %in% all_tips]
    
    if (length(tips) == 0) return(0)
    
    abund <- table(tips)
    spp   <- names(abund)
    
    if (length(spp) == 1) return(as.numeric(abund))
    
    sub_tree <- ape::keep.tip(phylo_tree, spp)
    
    edge_lengths <- sub_tree$edge.length
    weights <- rep(1, length(edge_lengths))
    
    for (k in seq_along(sub_tree$tip.label)) {
      sp <- sub_tree$tip.label[k]
      edge_idx <- which(sub_tree$edge[,2] == k)
      weights[edge_idx] <- abund[sp]
    }
    
    sum(edge_lengths * weights)
    
  },
  numeric(1)
)
beep()
# GRID
col_name <- paste0("AfaithPD_", window_radius_m, "m")
grid_sf[[col_name]] <- pd_vec

grid_vect  <- vect(grid_sf)
r_template <- rast(
  ext(grid_vect),
  resolution = 1,
  crs        = crs(grid_vect)
)
pd_rast <- rasterize(grid_vect, r_template, field = col_name)
plot(pd_rast)
out_raster <- file.path(output_folder,
                        paste0("AfaithPD_", window_radius_m, "m.tif"))
writeRaster(pd_rast, filename = out_raster, overwrite = TRUE)
#########################################################################
# MAP
fill_limits <- range(grid_sf[[col_name]], na.rm = TRUE)

pd_map <- ggplot() +
  geom_sf(
    data  = grid_sf,
    aes(fill = .data[[col_name]]),
    color = NA
  ) +
  scale_fill_gradientn(
    colours = viridis::viridis(5, option = "C"),
    limits  = fill_limits,
    oob     = scales::squish,
    name    = paste0("Faith's abundance PD\n(", window_radius_m, " m neighbourhood)\nMyr"),
    guide   = guide_colourbar(barwidth = 0.5, barheight = 15)
  ) +
  coord_sf(crs = st_crs(grid_sf)) +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    axis.title      = element_blank()
  ) +
  ggtitle(paste0(
    "Faith's abundance Phylogenetic Diversity (",
    window_radius_m,
    " m Moving Window, 1 m Grid)"
  ))

print(pd_map)
beep()
