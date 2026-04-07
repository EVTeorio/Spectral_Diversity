
library(terra)
library(sf)
library(dplyr)
library(snow)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

hsi_dir  <- "Spectral_Diversity/VegIndex_NA_trimmed"
out_dir  <- "Spectral_Diversity/Quad_Spectra/20m_VegIndex"
shp_path <- "Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp"

# Parallel
n_cores <- max(1, parallel::detectCores() - 2)

# Read in quads
quads <- st_read(shp_path, quiet = TRUE)
quads_sf <- quads %>% select(-matches("^Dscrptn"))
# quads_sf <- quads_sf %>%
#   mutate(Name = sub_id) %>%
#   select(-sub_id)    

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


# Starting place, If loop had to be interupted
#---------------------------------------
start_idx <- which(quads_sf$Name == "0")
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
  quads_sf <- st_read("Spectral_Diversity/Quad_Scale_SHPs/PR_20m.shp", quiet = TRUE)
  quads_sf <- quads_sf[, !grepl("^Dscrptn", names(quads_sf))]
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

#---------------------------------------------------------------------------------
# Run in parallel

results <- parLapply(cl, quad_indices, process_quad)
beep(3) #lets me know when its done

# Stop cluster
stopCluster(cl)


