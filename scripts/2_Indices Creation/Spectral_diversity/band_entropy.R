
library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

spec_dir <- "Quad_Spectra/20m_resampled_5nm"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Indices_SHPs/Spectral_diversitySHPs/SumBandEntropy_20m_masked_5nm.shp"

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

# --- Parameters ---
n_cores  <- max(1, parallel::detectCores() - 3)
best_band <- "603 nm"
threshold <- 0.01898618
direction <- "<"
n_bins <- 100

# --- Shannon entropy function ---
shannon_entropy <- function(x, n_bins) {
  h <- hist(x, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

# --- Per-band entropy function ---
per_band_entropy <- function(f) {
  
  r <- rast(f)
  
  # --- Apply masking logic ---
  band_mask <- r[[best_band]]
  sunlit_mask <- if (direction == ">") {
    band_mask < threshold
  } else {
    band_mask > threshold
  }
  
  r_masked <- mask(r, sunlit_mask, maskvalues = 0)
  
  # Extract matrix
  X <- values(r_masked, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) < 10) {
    return(tibble(
      Name = str_extract(basename(f), "\\d+"),
      spectral_entropy = NA_real_
    ))
  }
  
  # --- Compute entropy per band ---
  band_entropies <- apply(X, 2, shannon_entropy, n_bins = n_bins)
  
  # Mean entropy across bands
  mean_entropy <- mean(band_entropies, na.rm = TRUE)
  
  tibble(
    Name = str_extract(basename(f), "\\d+"),
    spectral_entropy = mean_entropy
  )
}

# --- Parallel processing ---
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, {library(terra); library(dplyr); library(stringr)})
clusterExport(cl, list("best_band", "threshold", "direction",
                       "n_bins", "shannon_entropy", "per_band_entropy"))

het_list <- parLapply(cl, ras_files, per_band_entropy)
beep(3)
stopCluster(cl)

het_df <- bind_rows(het_list)

# --- Merge with shapefile ---
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name) |>
  left_join(het_df, by = "Name")

plot(quads_sf)

st_write(quads_sf, out_shp, delete_layer = TRUE)
