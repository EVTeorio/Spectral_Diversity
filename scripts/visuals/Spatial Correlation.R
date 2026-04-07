

library(sf)
library(dplyr)
library(ggplot2)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

spectral_var <- st_read("Spectral_Diversity/Analysis_SHPs/PR_20m_spec_PCA_VAR.shp")

div_index <- st_read("Spectral_Diversity/Analysis_SHPs/sp_diversity_20m.shp")

head(div_index)
head(spectral_var)

#--------------------------------------------------
# 1. Join shapefiles by Name
#--------------------------------------------------

# Drop geometry from one layer to avoid conflicts
data_joined <- div_index %>%
  st_drop_geometry() %>%
  select(Name, shannon) %>%
  left_join(
    spectral_var %>%
      st_drop_geometry() %>%
      select(Name, spctrl_),
    by = "Name"
  )

# Remove any rows with missing values
data_joined <- data_joined %>%
  filter(!is.na(shannon), !is.na(spctrl_))

#--------------------------------------------------
# 2. Correlation analysis
#--------------------------------------------------

cor_test <- cor.test(
  data_joined$spctrl_,
  data_joined$shannon,
  method = "pearson" 
)

print(cor_test)

#--------------------------------------------------
# 3. Linear regression
#--------------------------------------------------

lm_model <- lm(shannon ~ spctrl_, data = data_joined)

summary(lm_model)

#--------------------------------------------------
# 4. Plot: Spectral variation vs Shannon diversity
#--------------------------------------------------

ggplot(data_joined, aes(x = spctrl_, y = shannon)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Spectral Variation (spctrl_)",
    y = "Shannon Species Diversity",
    title = "Relationship between Spectral Variation and Shannon Diversity"
  ) +
  theme_minimal()
