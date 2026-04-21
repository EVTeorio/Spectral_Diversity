
library(sf)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

spectral_div <- st_read("Indices_SHPs/Spectral_diversitySHPs/SA_entropy_20m_masked_5nm.shp")
plot(spectral_div$spctrl_)
plot(spectral_div)

# 1. Make valid
spectral_div <- st_make_valid(spectral_div)

# 2. Project to a planar CRS (UTM zone for your area in Alabama)
spectral_div_proj <- st_transform(spectral_div, 32616)  # WGS84 / UTM zone 16N

# 3. Snap to fix tiny gaps
spectral_div_proj <- st_snap(spectral_div_proj, spectral_div_proj, tolerance = 0.01)

# 4. Compute neighbors
neighbors <- st_touches(spectral_div_proj)

# 5. Count neighbors
spectral_div$n_neighbors <- lengths(neighbors)

# 6. Combine edge + corner
edge_polygons <- spectral_div[spectral_div$n_neighbors <= 5, ]
plot(st_geometry(spectral_div), col = "lightgrey")
plot(st_geometry(edge_polygons), col = "red", add = TRUE)

edge_polygons$Name
