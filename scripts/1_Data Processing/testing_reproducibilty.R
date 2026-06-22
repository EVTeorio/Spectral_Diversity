setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

origin_dir <- "Spectral_Diversity/Quad_Spectra/50m"
test_dir  <- "Spectral_Diversity/Quad_Spectra/50m_test"



library(terra)

# Directories
test_path <- "Spectral_Diversity/Quad_Spectra/50m_test"
orig_path <- "Spectral_Diversity/Quad_Spectra/50m"

# Get only the ENVI binary files (exclude hdr and aux.xml)
test_dir <- list.files(test_path, full.names = TRUE)
test_files <- test_dir[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", test_dir)]

orig_dir <- list.files(orig_path, full.names = TRUE)
orig_files <- orig_dir[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", orig_dir)]

# Match by filename
common_names <- intersect(basename(test_files), basename(orig_files))

# Sample 4 rasters
set.seed(1234)
sample_names <- sample(common_names, min(2, length(common_names)))

compare_raster_pair <- function(fname) {
  
  test_file <- test_files[basename(test_files) == fname]
  orig_file <- orig_files[basename(orig_files) == fname]
  
  r_test <- rast(test_file)
  r_orig <- rast(orig_file)
  
  # Compare metadata
  same_dims   <- all(dim(r_test) == dim(r_orig))
  same_extent <- ext(r_test) == ext(r_orig)
  same_crs    <- crs(r_test) == crs(r_orig)
  
  # Compare all values
  vals_test <- values(r_test, mat = FALSE)
  vals_orig <- values(r_orig, mat = FALSE)
  
  #same_values <- identical(vals_test, vals_orig)
  
  same_values <- isTRUE(
    all.equal(
      vals_test,
      vals_orig,
      tolerance = 1e-8,
      check.attributes = FALSE
    )
  )
  
  data.frame(
    raster = fname,
    same_dims = same_dims,
    same_extent = same_extent,
    same_crs = same_crs,
    same_values = same_values,
    identical = same_dims & same_extent & same_crs & same_values
  )
}

results <- do.call(
  rbind,
  lapply(sample_names, compare_raster_pair)
)

print(results)

























