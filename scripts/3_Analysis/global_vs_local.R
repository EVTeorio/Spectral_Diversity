

library(terra)
library(snow)
library(sf)
library(dplyr)
library(stringr)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")

masked <- st_read("Spectral_Diversity/Analysis_SHPs/global_PCA_VAR_20m_masked_5nm.shp")
global <- st_read("Spectral_Diversity/Analysis_SHPs/global_PCA_VAR_20m.shp")
sp_diversity <- st_read("Spectral_Diversity/Analysis_SHPs/sp_diversity_20m.shp")

head(local)
head(global)


x <- masked$spctrl_
y <- global$spctrl_

#comparing between local and global spectral heterogeneity 
cor(x, y, method = "spearman", use = "complete.obs")

#comparing between species diversity and local spectral heterogeneity 
cor(sp_diversity$shannon, x, method = "spearman", use = "complete.obs")

#comparing between species diversity and global spectral heterogeneity 
cor(sp_diversity$shannon, y, method = "spearman", use = "complete.obs")




