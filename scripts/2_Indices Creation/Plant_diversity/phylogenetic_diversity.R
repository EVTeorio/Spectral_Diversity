library(sf)
library(dplyr)
library(picante)
library(ape)
library(V.PhyloMaker2)
library(beepr)

# --- SETTINGS ---
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

shapefile_path <- "Indices_SHPs/species_wieghted_matrix_20m.shp"
taxa_path      <- "51sp_taxanomy.csv"
output_shp     <- "Indices_SHPs/Diversity_SHPs/Phylogenetic_Diversity_20m.shp"
crs_proj       <- 26916

# --- READ DATA ---
quads <- st_read("Quad_Scale_SHPs/PR_20m.shp")
quads <- quads %>% dplyr::select(-matches("^Dscrptn"))
quads <- st_transform(quads, crs_proj)
quads_sf <- st_read(shapefile_path)
quads_sf <- st_transform(quads_sf, crs_proj)

taxa_table <- read.csv(taxa_path)
taxa_table <- taxa_table %>%
  filter(sp_code != "COOB2")

# List of species columns
species_cols <- quads_sf %>%
  dplyr::select(-Name) %>%
  st_drop_geometry() %>%
  colnames()

# --- BUILD PHYLOGENETIC TREE ---
phylo_input <- taxa_table %>%
  mutate(species = paste(genus, species, sep = " ")) %>%
  dplyr::select(species, genus, family) %>%
  filter(species != "Carya sp") %>%      # Drop unknown species
  distinct(species, .keep_all = TRUE)

phylo_tree <- phylo.maker(
  sp.list   = phylo_input,
  tree      = GBOTB.extended.TPL,
  nodes     = nodes.info.1.TPL,
  scenarios = "S3"
)$scenario.3
beep()

# Map sp_code to tip labels (matching phylo tree)
sp_to_tip <- setNames(
  paste(taxa_table$genus, taxa_table$species, sep = "_"),
  taxa_table$sp_code
)
all_tips <- phylo_tree$tip.label

# --- FUNCTION TO CALCULATE FAITH'S PD PER POLYGON ---
calculate_faithPD <- function(abundances, species_cols) {
  present_species <- species_cols[abundances > 0]
  tip_labels <- sp_to_tip[present_species]
  tip_labels <- tip_labels[tip_labels %in% all_tips]
  
  if(length(tip_labels) == 0) return(0)
  if(length(tip_labels) == 1) return(1)
  
  sub_tree <- ape::keep.tip(phylo_tree, tip_labels)
  sum(sub_tree$edge.length)
}

# --- FUNCTION TO CALCULATE RAO PD ---
calculate_raoPD <- function(abundances, species_cols) {
  present_species <- species_cols[abundances > 0]
  tip_labels <- sp_to_tip[present_species]
  tip_labels <- tip_labels[tip_labels %in% all_tips]
  
  if(length(tip_labels) <= 1) return(0)
  
  abund <- table(tip_labels)
  sub_tree <- ape::keep.tip(phylo_tree, names(abund))
  dist_mat <- cophenetic(sub_tree)
  
  p <- as.numeric(abund[sub_tree$tip.label])
  p <- p / sum(p)
  
  sum(dist_mat * outer(p, p))
}

# --- FUNCTION TO CALCULATE ABUNDANCE-WEIGHTED FAITH PD ---
calculate_AfaithPD <- function(abundances, species_cols) {
  # Select species with positive abundance
  present_species <- species_cols[abundances > 0]
  tip_labels <- sp_to_tip[present_species]
  tip_labels <- tip_labels[tip_labels %in% all_tips]
  
  if(length(tip_labels) == 0) return(0)
  
  # Get actual abundances for the tips
  tip_abundances <- abundances[species_cols %in% present_species]
  names(tip_abundances) <- tip_labels
  
  if(length(tip_labels) == 1) return(as.numeric(tip_abundances))
  
  # Keep only relevant subtree
  sub_tree <- ape::keep.tip(phylo_tree, tip_labels)
  
  # Initialize edge weights (default = 1)
  edge_lengths <- sub_tree$edge.length
  weights <- rep(1, length(edge_lengths))
  
  # Multiply tip edges by **actual abundance**
  for (k in seq_along(sub_tree$tip.label)) {
    sp <- sub_tree$tip.label[k]
    edge_idx <- which(sub_tree$edge[,2] == k)  # tip edge index
    weights[edge_idx] <- tip_abundances[sp]
  }
  
  # Abundance-weighted PD
  sum(edge_lengths * weights)
}
# --- APPLY TO EACH POLYGON ---
abund_matrix <- quads_sf %>% st_drop_geometry() %>% dplyr::select(all_of(species_cols))

quads$faithPD  <- apply(abund_matrix, 1, calculate_faithPD, species_cols = species_cols)
quads$raoPD    <- apply(abund_matrix, 1, calculate_raoPD, species_cols = species_cols)
quads$AfaithPD <- apply(abund_matrix, 1, calculate_AfaithPD, species_cols = species_cols)

st_write(quads, output_shp)
