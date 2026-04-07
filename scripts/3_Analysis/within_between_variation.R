
library(terra)
library(snow)
library(dplyr)
library(stringr)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

spec_dir <- "Spectral_Diversity/Quad_Spectra/50m_resampled_5nm"
ras_files <- list.files(spec_dir, full.names = TRUE)
ras_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

# --- Parameters ---
n_cores <- max(1, parallel::detectCores() - 2)
best_band <- "603 nm"
threshold <- 0.01898618
direction <- "<"

# --- Worker function ---
tile_stats <- function(f) {
  
  r <- rast(f)
  
  # --- Apply masking logic ---
  band_mask <- r[[best_band]]
  sunlit_mask <- if (direction == ">") {
    band_mask < threshold
  } else {
    band_mask > threshold
  }
  
  r_masked <- mask(r, sunlit_mask, maskvalues = 0)
  
  X <- values(r_masked, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) < 10) {
    return(NULL)
  }
  
  list(
    name = str_extract(basename(f), "\\d+"),
    mean_spec = colMeans(X),
    within_var = mean(apply(X, 2, var))
  )
}

# --- Parallel ---
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, {library(terra); library(stringr)})
clusterExport(cl, list("best_band", "threshold", "direction", "tile_stats"))

tile_list <- parLapply(cl, ras_files, tile_stats)
beep(3)
stopCluster(cl)

# Remove NULL tiles
tile_list <- tile_list[!sapply(tile_list, is.null)]

# --- Extract results ---
tile_means <- do.call(rbind, lapply(tile_list, function(x) x$mean_spec))
within_vars <- sapply(tile_list, function(x) x$within_var)

# --- Compute variance components ---

# Mean within-tile variance
mean_within_variance <- mean(within_vars)

# Between-tile variance (variance of tile mean spectra)
between_variance <- mean(apply(tile_means, 2, var))

# Ratio
variance_ratio <- between_variance / mean_within_variance

# Output summary
cat("\nMean Within-Tile Variance:", mean_within_variance, "\n")
cat("Between-Tile Variance:", between_variance, "\n")
cat("Between / Within Ratio:", variance_ratio, "\n")
