

library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

spec_dir <- "Quad_Spectra/20m_resampled_5nm"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Indices_SHPs/Spectral_diversitySHPs/global_PCA_20m_masked_5nm.shp"

ras_files <- list.files(spec_dir, full.names = TRUE)
ras_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]
# --- Exclusion list ---
nums <- c(
  1424,1423,1422,1420,1421,1419,1418,1414,
  1521,1522,1523,1524,1520,1519,
  1624,1622,1623,1621,1620,
  1724,1723,1722,1721,
  1824,1823,1822,
  1923,1924,1922,1921,
  1322,1321,1319,1320,1318,1317,1316,1315,1314,1313,
  1221,1220,1219,1216,1215,1213,1214,1212,1211,
  1120,1119,1115,1114,1113,1112,1111,1110,
  1014,1013,1010,1009,
  909,908,24
)
# Filter out excluded quads
ras_files <- ras_files[!str_extract(basename(ras_files), "\\d+") %in% as.character(nums)]


# Parallel global PCA 
n_cores  <- max(1, parallel::detectCores() - 2)
n_sample <- 2000

cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, library(terra))

# EXPORT variables needed by workers
clusterExport(cl, list = "n_sample")

X_list <- parLapply(cl, ras_files, function(f) {
  
  r <- rast(f)
  
  X <- values(r, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) > n_sample) {
    X <- X[sample(nrow(X), n_sample), , drop = FALSE]
  }
  
  X
})
beep()
stopCluster(cl)

X_global <- do.call(rbind, X_list)

X_global_scaled <- scale(X_global)

global_center <- attr(X_global_scaled, "scaled:center")
global_scale  <- attr(X_global_scaled, "scaled:scale")

global_pca <- prcomp(X_global_scaled, center = FALSE, scale. = FALSE)
beep(3)

spectral_heterogeneity_global <- function(r, pca, center, scale,
                                          best_band = "603 nm",
                                          threshold = 0.01898618,
                                          direction = "<",
                                          n_pcs = 3) {
  
  # Select band for masking
  band_603 <- r[[best_band]]
  
  # Apply sunlit condition
  sunlit_mask <- if (direction == ">") {
    band_603 < threshold
  } else {
    band_603 > threshold
  }
  
  # Mask raster stack (keep only sunlit pixels)
  r_masked <- mask(r, sunlit_mask, maskvalues = 0)
  
  # Extract values
  X <- values(r_masked, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) < 10) return(NA_real_)
  
  # Apply global scaling
  Xs <- sweep(X, 2, center, "-")
  Xs <- sweep(Xs, 2, scale, "/")
  
  # Project into global PCA space
  scores <- Xs %*% pca$rotation[, 1:n_pcs, drop = FALSE]
  
  ctr <- colMeans(scores)
  mean(sqrt(rowSums((scores - ctr)^2)))
}

het_df <- bind_rows(lapply(ras_files, function(f) {
  
  message("Processing: ", basename(f))
  
  r <- rast(f)
  
  tibble(
    Name = str_extract(basename(f), "\\d+"),
    spectral_heterogeneity =
      spectral_heterogeneity_global(
        r,
        pca   = global_pca,
        center = global_center,
        scale  = global_scale,
        n_pcs  = 3
      )
  )
}))
beep(3)
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name) |>
  left_join(het_df, by = "Name")
plot(quads_sf)

st_write(quads_sf, out_shp, delete_layer = TRUE)











