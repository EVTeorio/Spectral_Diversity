
# Load library
library(terra)
library(sf)
library(dplyr)
library(lwgeom)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/2025_Biomass_LiDAR/")

# Read the KML
kml_sf <- st_read("Quad_Scale_SHPs/20ha.kml")

kml_sf <- st_zm(kml_sf, drop = TRUE, what = "ZM")


out_dir  <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/2025_Biomass_LiDAR/Quad_Scale_SHPs"
out_name <- "PR_20m"   

st_write(
  kml_sf,
  dsn = out_dir,
  layer = out_name,
  driver = "ESRI Shapefile")
  
#############################################3
#      Subsectin quads into 4 equal parts   #########
###     = 10m quads     ##############################

# Read the KML
sf_20m <- st_read("Quad_Scale_SHPs/PR_20m.shp")
head(sf_20m)

# ---- 2. Project to meters ----
sf_utm <- st_transform(sf_20m, 26916)

# ---- 3. Function to split a rotated square into 4 equal squares ----
split_rotated_square <- function(poly, name) {
  
  poly <- st_make_valid(poly)
  
  # Centroid
  ctr <- st_centroid(poly)
  
  # Oriented bounding box
  obb <- st_minimum_rotated_rectangle(poly)
  
  # Rectangle coordinates
  coords <- st_coordinates(obb)[, 1:2]
  
  # Orientation vectors
  v1 <- coords[2, ] - coords[1, ]
  v2 <- coords[3, ] - coords[2, ]
  
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  # Long enough to cut polygon
  bb <- st_bbox(poly)
  L <- max(bb$xmax - bb$xmin, bb$ymax - bb$ymin) * 2
  
  # Split lines
  line1 <- st_linestring(rbind(
    st_coordinates(ctr) - v1 * L,
    st_coordinates(ctr) + v1 * L
  ))
  
  line2 <- st_linestring(rbind(
    st_coordinates(ctr) - v2 * L,
    st_coordinates(ctr) + v2 * L
  ))
  
  split_lines <- st_sfc(line1, line2, crs = st_crs(poly))
  
  # Split polygon → sfc
  pieces <- st_split(poly, split_lines) |>
    st_collection_extract("POLYGON")
  
  # Convert to sf explicitly
  pieces_sf <- st_as_sf(pieces)
  
  # Add labels
  pieces_sf$sub_id <- paste0(name, "_", letters[1:nrow(pieces_sf)])
  
  pieces_sf
}

# ---- 4. Apply to all polygons ----
quarters <- do.call(
  rbind,
  lapply(seq_len(nrow(sf_utm)), function(i) {
    split_rotated_square(sf_utm[i, ], sf_utm$Name[i])
  })
)

# ---- 5. (Optional) Transform back to WGS84 ----
quarters_wgs84 <- st_transform(quarters, 4326)
warnings()

plot(quarters_wgs84)
st_write(quarters_wgs84, dsn = file.path(out_dir, paste0("PR_10m.shp")))
