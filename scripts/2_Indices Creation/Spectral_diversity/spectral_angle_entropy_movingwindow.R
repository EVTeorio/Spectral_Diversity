
library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

# --- Paths ---
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")
spec_dir <- "Quad_Spectra/20m_resampled_5nm"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Indices_SHPs/Spectral_diversitySHPs/SA_entropyMW_20m_masked_5nm.shp"

ras_files <- list.files(spec_dir, full.names = TRUE)
ras_files <- ras_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", ras_files)]

# # --- Exclusion list ---
# nums <- c(
#   1424,1423,1422,1420,1421,1419,1418,1414,
#   1521,1522,1523,1524,1520,1519,
#   1624,1622,1623,1621,1620,
#   1724,1723,1722,1721,
#   1824,1823,1822,
#   1923,1924,1922,1921,
#   1322,1321,1319,1320,1318,1317,1316,1315,1314,1313,
#   1221,1220,1219,1216,1215,1213,1214,1212,1211,
#   1120,1119,1115,1114,1113,1112,1111,1110,
#   1014,1013,1010,1009,
#   909,908,24
# )
# ras_files <- ras_files[!str_extract(basename(ras_files), "\\d+") %in% as.character(nums)]

# --- Parameters ---
n_cores        <- max(1, parallel::detectCores() - 2)
n_sample       <- 1000
best_band      <- "603 nm"
threshold      <- 0.01898618
direction      <- "<"
n_bins         <- 40
window_radius  <- 2   # metres — drives moving window size within each quad

# --- Functions ---
shannon_entropy <- function(x, n_bins) {
  h     <- hist(x, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

spectral_angle_entropy_window <- function(X, n_sample, n_bins) {
  # Normalize
  norms  <- sqrt(rowSums(X^2))
  valid  <- norms > 0
  X      <- X[valid, , drop = FALSE]
  norms  <- norms[valid]
  if (nrow(X) < 3) return(NA_real_)
  
  X_norm <- X / norms
  
  # Subsample
  n <- nrow(X_norm)
  if (n > n_sample) X_norm <- X_norm[sample(n, n_sample), , drop = FALSE]
  n <- nrow(X_norm)
  
  # Compare each pixel to mean spectrum (O(n) instead of O(n^2))
  mean_spec <- colMeans(X_norm)
  mean_spec <- mean_spec / sqrt(sum(mean_spec^2))
  dots      <- X_norm %*% mean_spec
  angles    <- acos(pmin(pmax(as.numeric(dots), -1), 1))
  
  shannon_entropy(angles, n_bins)
}

spectral_entropy_windowed <- function(f) {
  r <- rast(f)
  
  # Sunlit mask
  band_mask   <- r[[best_band]]
  sunlit_mask <- if (direction == ">") band_mask < threshold else band_mask > threshold
  r_masked    <- mask(r, sunlit_mask, maskvalues = 0)
  
  # Project to metric CRS for accurate metre-based buffering
  r_m <- project(r_masked, "EPSG:26916")
  
  # Build grid of centroids at native resolution
  r_template <- rast(ext(r_m), resolution = res(r_m), crs = crs(r_m))
  xy         <- as.data.frame(xyFromCell(r_template, 1:ncell(r_template)))
  
  # Convert centroids to sf and buffer by window_radius
  pts_sf  <- st_as_sf(xy, coords = c("x", "y"), crs = 26916)
  wins_sf <- st_buffer(pts_sf, dist = window_radius)
  
  # Compute entropy for each window
  win_entropies <- sapply(seq_len(nrow(wins_sf)), function(i) {
    win   <- vect(wins_sf[i, ])
    r_win <- crop(r_m, win)
    r_win <- mask(r_win, win)
    
    X <- values(r_win, mat = TRUE)
    X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
    X <- X[rowSums(X) > 0, , drop = FALSE]
    
    if (nrow(X) < 3) return(NA_real_)
    spectral_angle_entropy_window(X, n_sample = n_sample, n_bins = n_bins)
  })
  
  # Mean entropy across all windows = quad heterogeneity
  mean_entropy <- mean(win_entropies, na.rm = TRUE)
  
  tibble(
    Name             = str_extract(basename(f), "\\d+"),
    spectral_entropy = mean_entropy
  )
}

# --- Parallel processing ---
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, {
  library(terra)
  library(sf)
  library(dplyr)
  library(stringr)
  terraOptions(memfrac = 0.3, progress = 0)
})
clusterExport(cl, list("n_sample", "best_band", "threshold", "direction",
                       "n_bins", "window_radius", "shannon_entropy",
                       "spectral_angle_entropy_window", "spectral_entropy_windowed"))

het_list <- parLapply(cl, ras_files, spectral_entropy_windowed)
beep(3)
stopCluster(cl)

het_df <- bind_rows(het_list)

# --- Merge with shapefile ---
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name) |>
  left_join(het_df, by = "Name")

plot(quads_sf)
st_write(quads_sf, out_shp, delete_layer = TRUE)




