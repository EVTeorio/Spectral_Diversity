
library(terra)
library(sf)
library(dplyr)
library(snow)
library(beepr)
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")
hsi_dir  <- "Spectral_Diversity/VegIndex_NA_trimmed"
out_dir  <- "Spectral_Diversity/Quad_Spectra/10m_VegIndex"
shp_path <- "Spectral_Diversity/Quad_Scale_SHPs/PR_10m.shp"
# Parallel
n_cores <- max(1, parallel::detectCores() - 2)
# Read in quads
quads <- st_read(shp_path, quiet = TRUE)
quads_sf <- quads %>% select(-matches("^Dscrptn"))
quads_sf <- quads_sf %>%
  mutate(Name = sub_id) %>%
  select(-sub_id)    
# List HSI files
allfiles <- list.files(hsi_dir, full.names = TRUE)
hsi_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]
# CRS
hsi_template_path <- hsi_files[1]
hsi_template <- rast(hsi_template_path)
# Preload raster info
#----------------------------------------------
hsi_info <- lapply(hsi_files, function(f) {
  r <- rast(f)
  e <- ext(r)
  list(
    file   = f,
    extent = e,
    center = c((e$xmin + e$xmax)/2, (e$ymin + e$ymax)/2)
  )
})
# Starting place, if loop had to be interrupted
#---------------------------------------
start_idx <- which(quads_sf$Name == "0_a")
quad_indices <- start_idx:nrow(quads_sf)
# SNOW cluster
#----------------------------------------
cl <- makeCluster(n_cores, type = "SOCK")
# Export variables
clusterExport(cl, c("quad_indices", "out_dir", "hsi_files"), envir = environment())
# Load packages on threads
clusterEvalQ(cl, {
  library(terra)
  library(sf)
  library(dplyr)
  terraOptions(memfrac = 0.3, progress = 0)
  quads_sf <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_10m.shp", quiet = TRUE)
  quads_sf <- quads_sf[, !grepl("^Dscrptn", names(quads_sf))]
  quads_sf <- quads_sf %>% # only for 10m
    mutate(Name = sub_id) %>%
    select(-sub_id)
  NULL
})
# Processing function
#-----------------------------
process_quad <- function(i) {
  quad_sf <- quads_sf[i, ]
  quad_id <- quad_sf$Name
  
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
    x_center <- (e$xmin + e$xmax) / 2
    y_center <- (e$ymin + e$ymax) / 2
    
    # Skip if does not intersect
    if (!relate(quad_vect, as.polygons(e, crs=crs(r)), "intersects")) next
    
    r_quad <- crop(r, quad_vect)
    r_quad <- mask(r_quad, quad_vect)
    r_quad <- resample(r_quad, quad_template, method="near")
    
    if (all(is.na(values(r_quad[[1]])))) next
    
    # ------------------------------------
    # 2D EUCLIDEAN DISTANCE FROM RASTER CENTER
    # ------------------------------------
    x_raster <- init(quad_template, "x")
    y_raster <- init(quad_template, "y")
    dist_r <- sqrt((x_raster - x_center)^2 + (y_raster - y_center)^2)
    # Set no-data pixels to Inf so they never win which.min (use for raw Spectra)
    no_data <- app(r_quad, function(x) ifelse(all(is.na(x)), 1, NA))
    ## #####--Use for VEG INDEX ############################
    #no_data <- app(r_quad, function(x) ifelse(sum(is.na(x)) >= 30, 1, NA))
    ##########################################################################
    
    dist_r[no_data == 1] <- Inf
    # ------------------------------------
    
    spectral_stack[[length(spectral_stack)+1]] <- r_quad
    distance_stack[[length(distance_stack)+1]] <- dist_r
  }
  
  if (length(spectral_stack)==0) return(NULL)
  
  distances <- rast(distance_stack)
  winner <- app(distances, which.min)
  
  # Direct pixel-wise selection by winner index
  r_out <- spectral_stack[[1]]
  vals_out <- matrix(NA, nrow=ncell(r_out), ncol=nlyr(r_out))
  
  for (k in seq_along(spectral_stack)) {
    idx <- which(values(winner) == k)
    if (length(idx) == 0) next
    vals_out[idx, ] <- values(spectral_stack[[k]])[idx, ]
  }
  
  values(r_out) <- vals_out
  
  writeRaster(r_out, file.path(out_dir, quad_id), filetype="ENVI", overwrite=TRUE)
  gc()
  return(quad_id)
}
#---------------------------------------------------------------------------------
# Run in parallel
results <- parLapply(cl, quad_indices, process_quad)
beep(3)
# Stop cluster
stopCluster(cl)
