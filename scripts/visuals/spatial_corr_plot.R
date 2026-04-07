

# Load libraries
library(sf)
library(dplyr)
library(ggplot2)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

spectral_var <- st_read("Spectral_Diversity/Analysis_SHPs/PR_20m_spec_PCA_VAR.shp")

div_index <- st_read("Spectral_Diversity/Analysis_SHPs/sp_diversity_20m.shp")

#--------------------------------------------------
# 1. Join data by Name
#--------------------------------------------------

data_map <- div_index %>%
  left_join(
    spectral_var %>%
      st_drop_geometry() %>% 
      select(Name, spctrl_),
    by = "Name"
  )

#--------------------------------------------------
# 2. Run linear regression and calculate residuals
#--------------------------------------------------
lm_model <- lm(shannon ~ spctrl_, data = data_map)
data_map$residuals <- residuals(lm_model)

#--------------------------------------------------
# 3. Plot spatial maps
#--------------------------------------------------

# Shannon diversity
ggplot(data_map) +
  geom_sf(aes(fill = shannon), color = NA) +
  scale_fill_viridis_c(option = "viridis") +
  labs(
    fill = "Shannon Diversity",
    title = "Spatial Distribution of Shannon Species Diversity"
  ) +
  theme_minimal()

# Spectral variation
ggplot(data_map) +
  geom_sf(aes(fill = spctrl_), color = NA) +
  scale_fill_viridis_c(option = "magma") +
  labs(
    fill = "Spectral Variation",
    title = "Spatial Distribution of Spectral Variation"
  ) +
  theme_minimal()

# Residuals from regression
ggplot(data_map) +
  geom_sf(aes(fill = residuals), color = NA) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    fill = "Residuals",
    title = "Residuals: Shannon Diversity ~ Spectral Variation"
  ) +
  theme_minimal()
