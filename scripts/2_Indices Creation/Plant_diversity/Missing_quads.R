

library(sf)
library(dplyr)
library(stringr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

# --- Paths ---
in_dir  <- "Indices_SHPs/Diversity_SHPs"
out_dir <- "Indices_SHPs/sp_exclude_quads"

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

exclude_names <- as.character(nums)

# --- Process each shapefile ---
shp_files <- list.files(in_dir, pattern = "\\.shp$", full.names = TRUE)

cat("Found", length(shp_files), "shapefiles to process\n\n")

for (f in shp_files) {
  fname <- basename(f)
  cat("Processing:", fname, "\n")
  
  shp <- st_read(f, quiet = TRUE)
  
  n_before <- nrow(shp)
  shp_filtered <- shp %>%
    filter(!as.character(Name) %in% exclude_names)
  n_after  <- nrow(shp_filtered)
  
  cat("  Removed:", n_before - n_after, "polygons |",
      "Remaining:", n_after, "\n")
  
  out_path <- file.path(out_dir, fname)
  st_write(shp_filtered, out_path, delete_layer = TRUE, quiet = TRUE)
  cat("  Written to:", out_path, "\n\n")
}

cat("Done. All shapefiles written to:", out_dir, "\n")



