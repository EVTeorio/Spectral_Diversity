
library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

# --- Paths ---
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")
spec_dir <- "Quad_Spectra/20m_smoothed_5nm"
shp_path <- "Quad_Scale_SHPs/PR_20m.shp"
out_shp  <- "Indices_SHPs/Spectral_diversitySHPs/QuadInclSA_entropy_20m_smoothed_masked_5nm.shp"

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
n_cores <- max(1, parallel::detectCores() - 3)
n_sample <- 2000
best_band <- "563 nm"
threshold <- 0.0305476
direction <- "<"

# --- Functions ---
spectral_angle <- function(X) {
  # Normalize spectra
  X_norm <- X / sqrt(rowSums(X^2))
  n <- nrow(X_norm)
  
  # Subsample if too many pixels
  if (n > n_sample) {
    set.seed(42)
    X_norm <- X_norm[sample(n, n_sample), , drop = FALSE]
  }
  n <- nrow(X_norm)
  
  # Pairwise angles
  angles <- c()
  for (i in 1:(n-1)) {
    dots <- X_norm[i, ] %*% t(X_norm[(i+1):n, ])
    angles <- c(angles, acos(pmin(pmax(dots, -1), 1)))
  }
  angles
}

shannon_entropy <- function(x, n_bins = 40) {
  h <- hist(x, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)   # probabilities sum to 1
  probs <- probs[probs > 0]           # remove zeros
  -sum(probs * log(probs))
}

spectral_entropy_global <- function(f) {
  r <- rast(f) # formerly r
  
  # Uncomment to apply shadow mask
  # Mask based on sunlit condition
  band_mask <- r[[best_band]]
  sunlit_mask <- if (direction == ">") {
    band_mask < threshold
  } else {
    band_mask > threshold
  }
  r_masked <- mask(r, sunlit_mask, maskvalues = 0)

  # Extract values
  X <- values(r_masked, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  if (nrow(X) < 10) return(tibble(Name = str_extract(basename(f), "\\d+"),
                                  spectral_entropy = NA_real_))
  
  angles <- spectral_angle(X)
  entropy <- shannon_entropy(angles)
  
  tibble(Name = str_extract(basename(f), "\\d+"),
         spectral_entropy = entropy)
}

# --- Parallel processing ---
cl <- makeCluster(n_cores, type = "SOCK")
clusterEvalQ(cl, {library(terra); library(dplyr); library(stringr)})
clusterExport(cl, list("ras_files", "n_sample", "best_band", "threshold", "direction",
                       "spectral_angle", "shannon_entropy", "spectral_entropy_global"))

het_list <- parLapply(cl, ras_files, spectral_entropy_global)
beep()
stopCluster(cl)

het_df <- bind_rows(het_list)
write.csv(het_df, "Indices_SHPs/20m_SA_smooth_masked_7_11.csv")

# --- Merge with shapefile ---
quads_sf <- st_read(shp_path, quiet = TRUE) |>
  select(Name) |>
  left_join(het_df, by = "Name")

plot(quads_sf)
st_write(quads_sf, out_shp, delete_layer = TRUE)



