library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)
library(tidyr)

# --- Paths ---
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")
spec_dir <- "Quad_Spectra/20m_smoothed_5nm"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Indices_SHPs/Spectral_diversitySHPs/QuadInclSA_entropy_20m_smoothed_masked_5nm_boot.shp"
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
ras_files <- ras_files[!str_extract(basename(ras_files), "\\d+") %in% as.character(nums)]

# --- Parameters ---
n_cores    <- max(1, parallel::detectCores() - 3)
n_sample   <- 5000
n_boot     <- 70         # number of bootstrap iterations per quad
best_band  <- "563 nm"
threshold  <- 0.0305476
direction  <- "<"

# --- Functions ---
spectral_angle <- function(X_norm) {
  # Expects already-normalised matrix; no further subsampling here
  n <- nrow(X_norm)
  angles <- c()
  for (i in 1:(n - 1)) {
    dots   <- X_norm[i, ] %*% t(X_norm[(i + 1):n, , drop = FALSE])
    angles <- c(angles, acos(pmin(pmax(dots, -1), 1)))
  }
  angles
}

shannon_entropy <- function(x, n_bins = 40) {
  h     <- hist(x, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

spectral_entropy_boot <- function(f) {
  
  r <- rast(f)
  
  # --- Shadow mask ---
  band_mask  <- r[[best_band]]
  sunlit_mask <- if (direction == ">") band_mask < threshold else band_mask > threshold
  r_masked   <- mask(r, sunlit_mask, maskvalues = 0)
  
  # --- Extract clean pixels ---
  X     <- values(r_masked, mat = TRUE)
  X     <- X[rowSums(is.na(X)) == 0, ]
  X     <- X[rowSums(X) > 0, ]
  n_pix <- nrow(X)                        # total available pixels
  
  quad_name <- str_extract(basename(f), "\\d+")
  
  if (n_pix < 10) {
    # Return a tidy tibble with NAs for all bootstrap iterations
    out <- tibble(
      Name      = quad_name,
      n_pixels  = n_pix,
      boot_iter = seq_len(n_boot),
      spectral_entropy = NA_real_
    )
    return(out)
  }
  
  # Normalise once; resample inside loop
  X_norm <- X / sqrt(rowSums(X^2))
  
  # --- Bootstrap loop (no set.seed → different sample each iteration) ---
  entropies <- vapply(seq_len(n_boot), function(b) {
    idx    <- sample(n_pix, min(n_sample, n_pix), replace = FALSE)
    angles  <- spectral_angle(X_norm[idx, , drop = FALSE])
    shannon_entropy(angles)
  }, numeric(1))
  
  tibble(
    Name             = quad_name,
    n_pixels         = n_pix,
    boot_iter        = seq_len(n_boot),
    spectral_entropy = entropies
  )
}

# --- Parallel processing ---
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, { library(terra); library(dplyr); library(stringr) })
clusterExport(cl, list(
  "ras_files", "n_sample", "n_boot",
  "best_band", "threshold", "direction",
  "spectral_angle", "shannon_entropy", "spectral_entropy_boot"
))

boot_list <- parLapply(cl, ras_files, spectral_entropy_boot)
beep(3)
stopCluster(cl)

# --- Long-format dataframe  ---
boot_df_long <- bind_rows(boot_list)

# --- Wide-format dataframe (one row per quad, columns boot_1 … boot_10) ---
boot_df_wide <- boot_df_long |>
  pivot_wider(
    id_cols     = c(Name, n_pixels),
    names_from  = boot_iter,
    names_prefix = "boot_",
    values_from = spectral_entropy
  ) |>
  # Summary stats across the 10 iterations, computed row-wise
  rowwise() |>
  mutate(
    boot_mean   = mean(c_across(starts_with("boot_")),   na.rm = TRUE),
    boot_sd     = sd(c_across(starts_with("boot_")),     na.rm = TRUE),
    boot_cv     = boot_sd / boot_mean,
    boot_median = median(c_across(starts_with("boot_")), na.rm = TRUE),
    boot_min    = min(c_across(starts_with("boot_")),    na.rm = TRUE),
    boot_max    = max(c_across(starts_with("boot_")),    na.rm = TRUE)
  ) |>
  ungroup()

# --- Save both formats ---
write.csv(boot_df_long, "Indices_SHPs/20m_SA_smooth_masked_boot_long.csv", row.names = FALSE)
write.csv(boot_df_wide, "Indices_SHPs/20m_SA_entrop_boot100_results.csv", row.names = FALSE)

# --- Merge MEAN entropy with shapefile ---
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name) |>
  left_join(
    boot_df_wide |> select(Name, n_pixels, boot_mean, boot_sd, boot_cv),
    by = "Name"
  )

plot(quads_sf["boot_mean"])
st_write(quads_sf, out_shp, delete_layer = TRUE)