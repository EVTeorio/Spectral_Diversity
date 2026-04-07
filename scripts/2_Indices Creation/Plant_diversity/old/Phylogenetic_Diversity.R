
library(sf)
library(dplyr)
library(tidyr)
library(ape)
library(picante)
library(tibble)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
justforcrs <- rast("E:/Updated LiDAR/PRFPD_CHM_leafOff.tiff")
sp_table <- read.csv("SOFOR/big_data/heirchy table prfdp.csv")
#sp_table <- sp_table[,-1]
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
census <- read.csv("PRstem_DL.csv")
unique(census$sp)
# Filter which tree Size and or CPI
census <- filter(census, census$DBH.2024 >= 200)
#census <- filter(census, census$crown.position >= 200)

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

census_joined

# ---- 1. Clean census data ----
# Keep only grouping variable and species code
census_clean <- census_joined |>
  st_drop_geometry() |>
  select(Name, sp) |>
  filter(!is.na(sp))

# ---- 2. Clean taxonomy table ----
# Use sp_code as terminal taxa so it matches census data
sp_table_clean <- sp_table |>
  select(Class, Order, Family, Genus, sp_code) |>
  distinct() |>
  mutate(
    Class   = factor(Class),
    Order   = factor(Order),
    Family  = factor(Family),
    Genus   = factor(Genus),
    sp_code = factor(sp_code)
  )

# ---- 3. Build taxonomy-based phylogenetic tree ----
tax_tree <- as.phylo(
  ~ Class / Order / Family / Genus / sp_code,
  data = sp_table_clean
)

# Assign branch lengths (Grafen method)
tax_tree <- compute.brlen(tax_tree, method = "Grafen")

# ---- 3b. Root the tree and make binary ----
# Root at first tip and resolve multifurcations
tax_tree <- root(tax_tree, outgroup = tax_tree$tip.label[1], resolve.root = TRUE)
tax_tree <- multi2di(tax_tree)

# Check tree is rooted and binary
stopifnot(is.rooted(tax_tree))
stopifnot(is.binary.tree(tax_tree))

# ---- 4. Build community matrix (presence/absence) ----
comm_matrix <- xtabs(
  presence ~ Name + sp,
  data = transform(census_clean, presence = 1)
) |>
  as.matrix()

# Force presence/absence
comm_matrix[comm_matrix > 1] <- 1

# ---- 5. Match species between tree and community matrix ----
shared_species <- intersect(colnames(comm_matrix), tax_tree$tip.label)

# Keep only shared species
comm_matrix <- comm_matrix[, shared_species, drop = FALSE]
tax_tree <- drop.tip(tax_tree, setdiff(tax_tree$tip.label, shared_species))

# ---- 6. Calculate Faith's Phylogenetic Diversity ----
pd_results <- pd(
  samp = comm_matrix,
  tree = tax_tree,
  include.root = TRUE
)

# ---- 7. Format output ----
pd_results <- pd_results |>
  rownames_to_column("Name")


# ---- 1. Keep only Name and geometry ----
quads_simple <- quads |>
  select(Name, geometry)

# ---- 2. Join PD/FD metrics by Name ----
# Assuming your PD table is `pd_results` with columns: Name, PD, SR, etc.
quads_fd <- quads_simple |>
  left_join(pd_results, by = "Name")

st_write(quads)









####################################################################3

#############################
library(ape)
library(phytools)
library(picante)

# ======================================
# Step 1: Explicit species list
# ======================================
species_cols <- intersect(
  as.character(heirachy$sp_code),
  colnames(quads_species)
)

stopifnot(length(species_cols) > 1)

# ======================================
# Step 2: Build species matrix
# ======================================
species_matrix <- quads_species %>%
  st_drop_geometry() %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# Force character colnames (CRITICAL)
colnames(species_matrix) <- as.character(colnames(species_matrix))

# Remove empty quads
species_matrix <- species_matrix[rowSums(species_matrix) > 0, , drop = FALSE]

# ======================================
# Step 3: Prepare hierarchy table
# ======================================
heirachy_clean <- heirachy %>%
  filter(sp_code %in% species_cols) %>%
  mutate(
    Class   = factor(Class),
    Order   = factor(Order),
    Family  = factor(Family),
    Genus   = factor(Genus),
    sp_code = factor(sp_code)
  )

# ======================================
# Step 4: Build pseudo-phylogenetic tree
# ======================================
tax_formula <- ~ Class / Order / Family / Genus / sp_code
phylo_tree <- as.phylo(tax_formula, data = heirachy_clean)

# Force character tip labels (CRITICAL)
phylo_tree$tip.label <- as.character(phylo_tree$tip.label)

# ======================================
# Step 5: Root the tree
# ======================================
phylo_tree <- root(
  phylo_tree,
  outgroup = phylo_tree$tip.label[1],
  resolve.root = TRUE
)

# ======================================
# Step 6: Assign branch lengths
# ======================================
phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))

# ======================================
# Step 7: Explicitly prune tree to matrix species
# ======================================
common_species <- intersect(
  colnames(species_matrix),
  phylo_tree$tip.label
)

stopifnot(length(common_species) > 0)

phylo_tree_pruned <- keep.tip(phylo_tree, common_species)

# ======================================
# Step 8: Run Faith's PD
# ======================================
pd_values <- pd(
  samp = species_matrix[, common_species, drop = FALSE],
  tree = phylo_tree_pruned,
  include.root = TRUE
)

# ======================================
# Step 9: Attach PD back to sf object
# ======================================
quads_species$Faith_PD <- NA_real_

quads_species$Faith_PD[
  match(rownames(pd_values), rownames(quads_species))
] <- pd_values$PD



