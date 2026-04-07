
library(terra)
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(beepr)
library(stringr)
library(snow)

# SETTINGS
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
hsi_path        <- "Spectral_Diversity/Quad_Spectra/20m_resampled_5nm"
output_folder   <- "Spectral_Diversity/Indices_SHPs/Diversity_Rasters"
crs_proj        <- 26916
window_radius_m <- 20
n_sample        <- 500
best_band       <- "603 nm"
threshold       <- 0.01898618
direction       <- "<"
n_bins          <- 40
n_cores         <- max(1, parallel::detectCores() - 3)

# --- Functions ---
spectral_angle_entropy <- function(X, n_sample, n_bins) {
  norms  <- sqrt(rowSums(X^2))
  X_norm <- X / norms
  
  n <- nrow(X_norm)
  if (n > n_sample) X_norm <- X_norm[sample(n, n_sample), , drop = FALSE]
  n <- nrow(X_norm)
  
  if (n < 3) return(NA_real_)
  
  angles <- c()
  for (i in 1:(n - 1)) {
    dots   <- X_norm[i, ] %*% t(X_norm[(i + 1):n, , drop = FALSE])
    angles <- c(angles, acos(pmin(pmax(as.numeric(dots), -1), 1)))
  }
  
  h     <- hist(angles, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

process_window <- function(i) {
  win   <- vect(window_sf[i, ])
  r_win <- crop(hsi_masked, win)
  r_win <- mask(r_win, win)
  
  X <- values(r_win, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) < 3) return(NA_real_)
  
  spectral_angle_entropy(X, n_sample = n_sample, n_bins = n_bins)
}

# LOAD & MOSAIC HSI
hsi_files <- list.files(hsi_path, full.names = TRUE)
hsi_files <- hsi_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", hsi_files)]

cat("Mosaicking", length(hsi_files), "HSI files...\n")
hsi_list <- lapply(hsi_files, rast)
hsi_rast <- do.call(mosaic, c(hsi_list, fun = "mean"))
beep()

# REPROJECT
hsi_rast <- project(hsi_rast, paste0("EPSG:", crs_proj))

# SUNLIT MASK
band_mask   <- hsi_rast[[best_band]]
sunlit_mask <- if (direction == ">") band_mask < threshold else band_mask > threshold
hsi_masked  <- mask(hsi_rast, sunlit_mask, maskvalues = 0)

# Write masked raster to temp file so workers can read it
tmp_hsi <- tempfile(fileext = ".tif")
writeRaster(hsi_masked, tmp_hsi, overwrite = TRUE)
cat("Masked HSI written to temp file\n")

# BUILD 1 x 1 m GRID
bbox    <- ext(hsi_masked)
grid_1m <- st_make_grid(
  x        = st_as_sfc(st_bbox(c(xmin = bbox$xmin, xmax = bbox$xmax,
                                 ymin = bbox$ymin, ymax = bbox$ymax),
                               crs = st_crs(crs_proj))),
  cellsize = 1,
  square   = TRUE
)
grid_sf <- st_sf(geometry = grid_1m, crs = st_crs(crs_proj)) %>%
  mutate(cell_id = row_number())

cat("Grid built:", nrow(grid_sf), "cells\n")

# MOVING WINDOW
centroids <- st_centroid(grid_sf)
window_sf <- st_buffer(centroids, dist = window_radius_m)
win_indices <- seq_len(nrow(window_sf))
beep()

# PARALLEL CLUSTER
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, {
  library(terra)
  library(sf)
  library(dplyr)
  terraOptions(memfrac = 0.3, progress = 0)
})
clusterExport(cl, list("window_sf", "tmp_hsi", "n_sample", "n_bins",
                       "spectral_angle_entropy", "process_window"),
              envir = environment())

# Load masked raster on each worker from temp file
clusterEvalQ(cl, {
  hsi_masked <- rast(tmp_hsi)
})

cat("Computing spectral angle entropy in parallel...\n")
sa_list <- parLapply(cl, win_indices, process_window)
beep(3)
stopCluster(cl)

# Clean up temp file
file.remove(tmp_hsi)

sa_vec <- unlist(sa_list)

# RASTERIZE
col_name <- paste0("SA_entropy_", window_radius_m, "m")
grid_sf[[col_name]] <- sa_vec

grid_vect  <- vect(grid_sf)
r_template <- rast(ext(grid_vect), resolution = 1, crs = crs(grid_vect))
sa_rast    <- rasterize(grid_vect, r_template, field = col_name)

out_raster <- file.path(output_folder,
                        paste0("SA_entropy_", window_radius_m, "m.tif"))
writeRaster(sa_rast, filename = out_raster, overwrite = TRUE)
cat("Raster written to:", out_raster, "\n")

# MAP
fill_limits <- range(grid_sf[[col_name]], na.rm = TRUE)

sa_map <- ggplot() +
  geom_sf(
    data  = grid_sf,
    aes(fill = .data[[col_name]]),
    color = NA
  ) +
  scale_fill_gradientn(
    colours = viridis::viridis(5, option = "C"),
    limits  = fill_limits,
    oob     = scales::squish,
    name    = paste0("SA Entropy\n(", window_radius_m, " m neighbourhood)"),
    guide   = guide_colourbar(barwidth = 0.5, barheight = 15)
  ) +
  coord_sf(crs = st_crs(grid_sf)) +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    axis.title      = element_blank()
  ) +
  ggtitle(paste0("Spectral Angle Entropy (",
                 window_radius_m, " m Moving Window, 1 m Grid)"))

print(sa_map)
beep()

