library(sf)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

csv_data <- read.csv("Indices_SHPs/20m_spectral_sp.csv")
  #dplyr::select(-X)

quads_with_elevation <- st_read("Indices_SHPs/Other_variables/20m_elevation.shp")
quads_with_elevation<- as.data.frame(quads_with_elevation) %>%
  dplyr::select(-geometry)
quads_with_elevation$Name <- as.integer(quads_with_elevation$Name)

merged_data <- csv_data %>%
  left_join(quads_with_elevation %>% select(Name, avg_lvt, elvtn_r), by = "Name") %>%
  dplyr::select(everything(),avg_lvt, elvtn_r, -geometry, geometry)

# # Merge the CSV with the shapefile based on the "Name" column
# merged_data <- csv_data %>%
#   left_join(quads_with_elevation %>% dplyr::select(Name, dmnnt_v), by = "Name") %>%
#   dplyr::select(everything(), dmnnt_v, -geometry, geometry)

# Write the updated data back to a new CSV
write.csv(merged_data, "Indices_SHPs/20m_spectral_sp.csv", row.names = FALSE) 
head(merged_data)
