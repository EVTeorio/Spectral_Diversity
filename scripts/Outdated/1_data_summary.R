
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/2025_Biomass_LiDAR/SOFOR")


summary <- read.csv("big_data/PaintRock_Spectra_Summary.CSV")
unique(summary$SpeciesID)
unique(summary$TreeID)


tree_counts <- summary %>%
  group_by(SpeciesID) %>%
  summarise(UniqueTreeCount = n_distinct(TreeID))

sum(tree_counts$UniqueTreeCount)

tree_counts_genus <- summary %>%
  group_by(Genus, SpeciesID) %>%
  summarise(UniqueTreeCount = n_distinct(TreeID),
            .groups = "drop")

