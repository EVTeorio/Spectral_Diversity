
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

# Read data
census <- read.csv("Spectral_Diversity/PR_tree_DL.csv")
csv <- read.csv("Spectral_Diversity/PR_tree_DL.csv")

csv <- filter(csv, census$DBH.2024 >= 200)
census <- filter(census, DBH.2024 >= 200 | crown.position %in% c(4, 5))

extra_in_census <- anti_join(census, csv, by = "StemTag")

species_summary <- census %>%
  group_by(sp) %>%
  summarise(
    n_individuals = n_distinct(StemTag)
  ) %>%
  arrange(desc(n_individuals))

write.csv(species_summary, 
          "Spectral_Diversity/species_summary.csv", 
          row.names = FALSE)
########################################################3

tree_df
tree_df$sp[tree_df$sp == "CACA18"] <- "CACA38"

write.csv(tree_df)