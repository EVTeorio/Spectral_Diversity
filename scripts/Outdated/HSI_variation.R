

library(terra)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

# -----------------------------
# Paths
# -----------------------------
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

spec_dir <- "Spectral_Diversity/Quad_Spectra/20m"
shp_path <- "Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Spectral_Diversity/Analysis_SHPs/PR_20m_spec_PCA_VAR.shp"

# -----------------------------
# Read quadrats
# -----------------------------
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name)

# -----------------------------
# Spectral heterogeneity function
# -----------------------------
spectral_heterogeneity <- function(r, n_pcs = 10) {
  
  X <- terra::values(r, mat = TRUE)
  
  # Drop NA / empty pixels
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) < 10) return(NA_real_)
  
  Xs <- scale(X)
  
  pca <- prcomp(Xs, center = TRUE, scale. = FALSE)
  
  scores <- pca$x[, 1:n_pcs, drop = FALSE]
  ctr <- colMeans(scores)
  
  mean(sqrt(rowSums((scores - ctr)^2)))
}

# -----------------------------
# Loop over rasters
# -----------------------------
ras_files <- list.files(spec_dir, full.names = TRUE)
ras_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

het_df <- bind_rows(lapply(ras_files, function(f) {
  
  message("Processing: ", basename(f))
  
  r <- rast(f)
  
  tibble(
    Name = str_extract(basename(f), "\\d+"),
    spectral_heterogeneity = spectral_heterogeneity(r, n_pcs = 10)
  )
}))

beep(3)
# -----------------------------
# Join to polygons
# -----------------------------
quads_sf <- quads_sf |>
  left_join(het_df, by = "Name")

# -----------------------------
# Write output
# -----------------------------
st_write(quads_sf, out_shp, delete_layer = TRUE)





