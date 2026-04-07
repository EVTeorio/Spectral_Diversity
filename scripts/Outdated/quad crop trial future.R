
#==================================================
# SNOW-parallel processing of quads with HSI tiles
#==================================================

library(terra)
library(sf)
library(dplyr)
library(snow)
library(beepr)

#--------------------------------------------------
# Settings
#--------------------------------------------------
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

hsi_dir  <- "Spectral_Diversity/HSI_NA_trimmed"
out_dir  <- "Spectral_Diversity/Quad_Spectra/20m"
shp_path <- "Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp"

# Parallel
n_cores <- max(1, parallel::detectCores() - 2)

#--------------------------------------------------
# Read and clean quads
#--------------------------------------------------
quads <- st_read(shp_path, quiet = TRUE)
quads_sf <- quads %>% select(-matches("^Dscrptn"))
# quads_sf <- quads_sf %>%
#   mutate(Name = sub_id) %>%   # replace Name with sub_id
#   select(-sub_id)    

#--------------------------------------------------
# List HSI files (exclude sidecars)
#--------------------------------------------------
allfiles <- list.files(hsi_dir, full.names = TRUE)
hsi_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]

#--------------------------------------------------
# CRS / template
#--------------------------------------------------
hsi_template_path <- hsi_files[1]
hsi_template <- rast(hsi_template_path)

#--------------------------------------------------
# Preload lightweight raster info only
#--------------------------------------------------
hsi_info <- lapply(hsi_files, function(f) {
  r <- rast(f)
  e <- ext(r)
  list(
    file   = f,
    extent = e,
    center = c((e$xmin + e$xmax)/2, (e$ymin + e$ymax)/2)
  )
})

#--------------------------------------------------
# Quad indices (example: start at sub50_33)
#--------------------------------------------------
start_idx <- which(quads_sf$Name == "0")
if (length(start_idx) == 0) stop("Quad not found in shapefile")
quad_indices <- start_idx:nrow(quads_sf)

#--------------------------------------------------
# SNOW cluster
#--------------------------------------------------
#--------------------------------------------------
# Minimal export to workers
#--------------------------------------------------
cl <- makeCluster(n_cores, type = "SOCK")

# Export only necessary simple variables
clusterExport(cl, c("quad_indices", "out_dir", "hsi_files"), envir = environment())

# Load packages on workers
clusterEvalQ(cl, {
  library(terra)
  library(sf)
  library(dplyr)
  terraOptions(memfrac = 0.3, progress = 0)
  # Each worker reads quads itself
  quads_sf <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp", quiet = TRUE)
  quads_sf <- quads_sf[, !grepl("^Dscrptn", names(quads_sf))]
  NULL
})

#--------------------------------------------------
# Processing function (workers read rasters themselves)
#--------------------------------------------------
process_quad <- function(i) {
  quad_sf <- quads_sf[i, ]
  quad_id <- quad_sf$Name
  
  # Convert to vect
  quad_vect <- vect(quad_sf)
  
  # Template raster
  hsi_template <- rast(hsi_files[1])
  quad_vect <- project(quad_vect, hsi_template)
  quad_template <- rast(ext(quad_vect), resolution=res(hsi_template), crs=crs(hsi_template))
  
  # Loop through HSI files
  spectral_stack <- list()
  distance_stack <- list()
  
  for (f in hsi_files) {
    r <- rast(f)
    e <- ext(r)
    center <- c((e$xmin+e$xmax)/2, (e$ymin+e$ymax)/2)
    
    # Skip if does not intersect
    if (!relate(quad_vect, as.polygons(e, crs=crs(r)), "intersects")) next
    
    r_quad <- crop(r, quad_vect)
    r_quad <- mask(r_quad, quad_vect)
    r_quad <- resample(r_quad, quad_template, method="near")
    
    if (all(is.na(values(r_quad[[1]])))) next
    
    dist_r <- distance(quad_template, vect(matrix(center, ncol=2), type="points", crs=crs(quad_template)))
    
    spectral_stack[[length(spectral_stack)+1]] <- r_quad
    distance_stack[[length(distance_stack)+1]] <- dist_r
  }
  
  if (length(spectral_stack)==0) return(NULL)
  
  distances <- rast(distance_stack)
  winner <- app(distances, which.min)
  
  r_out <- spectral_stack[[1]]
  values(r_out) <- NA
  
  for (k in seq_along(spectral_stack)) {
    r_out <- cover(r_out, mask(spectral_stack[[k]], winner==k))
  }
  
  writeRaster(r_out, file.path(out_dir, quad_id), filetype="ENVI", overwrite=TRUE)
  gc()
  return(quad_id)
}

#--------------------------------------------------
# Run in parallel
#--------------------------------------------------
results <- parLapply(cl, quad_indices, process_quad)
beep()
# Stop cluster
stopCluster(cl)

##################
#######################################
######################################
############

#==================================================
# Parallel HSI Quad Extraction (SNOW, Windows-safe)
#==================================================

#------------------------------
# Libraries
#------------------------------
library(terra)
library(sf)
library(parallel)

#------------------------------
# Settings
#------------------------------
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

# Input/output directories
hsi_dir  <- "Spectral_Diversity/HSI_NA_trimmed"
out_dir  <- "Spectral_Diversity/Quad_Spectra/50m"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of cores (leave 2 free)
n_cores <- max(1, parallel::detectCores() - 2)
message("Using ", n_cores, " cores for parallel processing")

#------------------------------
# Read quads
#------------------------------
quads_sf <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_50m.shp", quiet = TRUE)
quads_sf <- quads_sf[, !grepl("^Dscrptn", names(quads_sf))]  # drop description fields
quad_indices <- seq_len(nrow(quads_sf))

#------------------------------
# List HSI rasters
#------------------------------
allfiles <- list.files(hsi_dir, full.names = TRUE)
hsi_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]
if (length(hsi_files) == 0) stop("No HSI files found in ", hsi_dir)

# Path to template raster for CRS/resolution
hsi_template_path <- hsi_files[1]

#------------------------------
# Setup SNOW cluster
#------------------------------
cl <- makeCluster(n_cores, type = "SOCK")

# Export small objects only
clusterExport(cl, c("quad_indices", "hsi_files", "out_dir", "hsi_template_path"), envir = environment())

# Load libraries and read quads on workers
clusterEvalQ(cl, {
  library(terra)
  library(sf)
  terraOptions(memfrac = 0.3, progress = 0)
  
  # Each worker reads quads locally
  quads_sf <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_50m.shp", quiet = TRUE)
  quads_sf <- quads_sf[, !grepl("^Dscrptn", names(quads_sf))]
  NULL
})

#------------------------------
# Quad processing function
#------------------------------
process_quad <- function(i) {
  
  quad_sf <- quads_sf[i, ]
  quad_id <- quad_sf$Name
  message("Processing quad: ", quad_id)
  
  # Convert to terra vector
  quad_vect <- vect(quad_sf)
  
  # Template raster for CRS and resolution
  hsi_template <- rast(hsi_template_path)
  quad_vect <- project(quad_vect, hsi_template)
  quad_template <- rast(ext(quad_vect), resolution=res(hsi_template), crs=crs(hsi_template))
  
  # Initialize stacks
  spectral_stack <- list()
  distance_stack <- list()
  
  # Loop through HSI tiles
  for (f in hsi_files) {
    r <- rast(f)
    e <- ext(r)
    center <- c((e$xmin + e$xmax)/2, (e$ymin + e$ymax)/2)
    
    # Skip if tile does not intersect quad
    if (!relate(quad_vect, as.polygons(e, crs=crs(r)), "intersects")) next
    
    # Crop, mask, resample
    r_quad <- crop(r, quad_vect)
    r_quad <- mask(r_quad, quad_vect)
    r_quad <- resample(r_quad, quad_template, method="near")
    
    if (all(is.na(values(r_quad[[1]])))) next  # skip if empty
    
    # Distance raster from tile center
    dist_r <- distance(quad_template, vect(matrix(center, ncol=2), type="points", crs=crs(quad_template)))
    
    spectral_stack[[length(spectral_stack)+1]] <- r_quad
    distance_stack[[length(distance_stack)+1]] <- dist_r
  }
  
  if (length(spectral_stack) == 0) return(NULL)  # skip if no data
  
  # Determine closest tile per pixel
  distances <- rast(distance_stack)
  winner <- app(distances, which.min)
  
  r_out <- spectral_stack[[1]]
  values(r_out) <- NA
  
  for (k in seq_along(spectral_stack)) {
    r_out <- cover(r_out, mask(spectral_stack[[k]], winner == k))
  }
  
  # Write output raster
  writeRaster(r_out, file.path(out_dir, quad_id), filetype="ENVI", overwrite=TRUE)
  gc()
  return(quad_id)
}

#------------------------------
# Run SNOW parallel
#------------------------------
results <- parLapply(cl, quad_indices, process_quad)

#------------------------------
# Stop cluster
#------------------------------
stopCluster(cl)

message("Processing complete. Quads processed: ", sum(!sapply(results, is.null)))




