install.packages("gsl")

library(dplyr)
library(ggplot2)
library(plotly)
library(vegan)


# Read data
summary <- read.csv("big_data/PaintRock_Spectra_Summary.CSV")

# Identify spectral columns (all numeric, excluding metadata)
spectra_cols <- summary %>%
  select(where(is.numeric)) %>%
  select(-X.1, -X, -ImageID, -TreeID) %>%
  colnames()

spectra <- summary[, spectra_cols]

# Bray–Curtis distance (or euclidean if scaled)
spectra_scaled <- scale(spectra)
dist_mat <- dist(spectra_scaled, method = "euclidean")

# PERMANOVA
adonis_res <- adonis2(
  dist_mat ~ SpeciesID,
  data = summary,
  permutations = 999,
  by = "margin"
)

# PERMANOVA
adonis_res <- adonis2(
  dist_mat ~ Genus,
  data = summary,
  permutations = 999,
  by = "margin"
)


adonis_res



####################################

pca <- prcomp(spectra, center = TRUE, scale. = TRUE)

summary(pca)


pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  TreeID = summary$TreeID,
  SpeciesID = summary$SpeciesID,
  Genus = summary$Genus
)


p_species <- plot_ly(
  data = pca_df,
  x = ~PC1,
  y = ~PC2,
  color = ~SpeciesID,
  type = "scatter",
  mode = "markers",
  marker = list(size = 7, opacity = 0.8),
  text = ~paste(
    "TreeID:", TreeID,
    "<br>SpeciesID:", SpeciesID,
    "<br>Genus:", Genus
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "Interactive PCA Colored by SpeciesID",
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2")
  )

p_species

#Genus
p_genus <- plot_ly(
  data = pca_df,
  x = ~PC1,
  y = ~PC2,
  color = ~Genus,
  type = "scatter",
  mode = "markers",
  marker = list(size = 7, opacity = 0.8),
  text = ~paste(
    "TreeID:", TreeID,
    "<br>SpeciesID:", SpeciesID
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "Interactive PCA Colored by Genus",
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2")
  )

p_genus

################################################

# Distance matrices
spec_dist <- dist(spectra_scaled, method = "euclidean")
species_dist <- dist(as.numeric(as.factor(summary$SpeciesID)))
genus_dist <- dist(as.numeric(as.factor(summary$Genus)))

# Mantel tests
mantel_species <- mantel(spec_dist, species_dist, permutations = 999)
mantel_genus <- mantel(spec_dist, genus_dist, permutations = 999)

mantel_species
mantel_genus
