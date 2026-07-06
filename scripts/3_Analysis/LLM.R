model_lm <- lm(spectral_het ~ pd * poly(elevation, 2), data = model_df)

# library(splines)
# 
# model_lm <- lm(spectral_het ~ pd * ns(elevation, df = 3), data = model_df)
# summary(model_lm)
# 
# model_df <- model_df %>%
#   mutate(
#     elev_bin = cut(
#       elevation,
#       breaks = quantile(elevation, probs = c(0, 0.33, 0.66, 1)),
#       labels = c("low", "mid", "high"),
#       include.lowest = TRUE
#     )
#   )
# 
# model_bin <- lm(spectral_het ~ pd * elev_bin, data = model_df)
# 
# summary(model_bin)


# =============================================================================
# 0. LIBRARIES
# =============================================================================
library(dplyr)
library(ggplot2)
library(sf)
library(janitor)
library(interactions)
library(car)        # VIF
library(lme4)       # LMM
library(lmerTest)   # p-values for lmer
library(performance) # optional diagnostics

# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

df <- read.csv("Quad_Values/20m_spectral_sp.csv") %>%
  clean_names() %>%
  filter(complete.cases(.))

df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)

# =============================================================================
# 2. VARIABLE DEFINITIONS
# =============================================================================

# Spectral heterogeneity (response)
entropy_var <- "sa_entropy_20m_masked_5nm"

# Elevation
elevation_var <- "avg_lvt"

# Phylogenetic diversity
pd_var <- "afaith_pd"

# Keep only needed columns
model_df <- df %>%
  select(all_of(c(entropy_var, elevation_var, pd_var))) %>%
  na.omit()

# Rename for convenience
model_df <- model_df %>%
  rename(
    spectral_het = all_of(entropy_var),
    elevation = all_of(elevation_var),
    pd = all_of(pd_var)
  )

# =============================================================================
# 3. CHECK COLLINEARITY
# =============================================================================

# Fit simple linear model for VIF check (no random effects here)
lm_check <- lm(spectral_het ~ pd + elevation, data = model_df)

vif_values <- vif(lm_check)
print(vif_values)

# Optional: correlation matrix
cor_matrix <- cor(model_df)
print(cor_matrix)

# =============================================================================
# 4. FIT LINEAR MIXED MODEL
# =============================================================================

# ---- IMPORTANT ----
# If you DO NOT have grouping structure, use lm() instead.
# If you DO, uncomment and modify random effects below.

# ---- OPTION A: No random effects (simple model) ----
model_lm <- lm(spectral_het ~ pd * elevation, data = model_df)
summary(model_lm)

# ---- OPTION B: With random effects (EDIT THIS) ----
# Example if you have a column like "site" or "plot":
model_df <- model_df %>%
  mutate(name = df$name) %>%
  select(name, everything())

model_lmm <- lmer(spectral_het ~ pd * elevation + (1 | name), data = model_df)

# summary(model_lmm)

# =============================================================================
# 5. MODEL DIAGNOSTICS
# =============================================================================

# Residual plots
par(mfrow = c(2,2))
plot(model_lm)

# Optional (better diagnostics)
check_model(model_lm)

# =============================================================================
# 6. VISUALIZATION
# =============================================================================

ggplot(model_df, aes(x = pd, y = spectral_het, color = elevation)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    x = "Phylogenetic Diversity",
    y = "Spectral Heterogeneity",
    color = "Elevation"
  )

# Interaction visualization
library(interactions)

interact_plot(model_lm, pred = pd, modx = elevation)
