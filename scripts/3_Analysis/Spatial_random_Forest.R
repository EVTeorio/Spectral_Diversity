
library(dplyr)
library(ranger)
library(adespatial)
library(sf)
library(ggplot2)
library(tibble)
library(spdep)


setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
df <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  filter(complete.cases(.))

# --- Setup ---
spectral_cols <- c("global_PCA_20m_masked_5nm", "SA_entropy_20m_masked_5nm",
                   "SumBandEntropy_20m_masked_5nm")

species_start <- which(names(df) == "CAGL8")
species_cols  <- names(df)[species_start:which(names(df) == "COOB2")]

# --- Convert geometry to sf and extract centroids ---
df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)

centroids <- df_sf %>%
  st_transform(26916) %>%
  st_centroid()

coords <- st_coordinates(centroids) %>%
  as.data.frame() %>%
  rename(x = X, y = Y)

# --- Build modeling dataframe ---
df_model <- df %>%
  select(Name, all_of(spectral_cols), all_of(species_cols)) %>%
  bind_cols(coords)

# --- Compute MEM spatial eigenvectors ---
cat("Computing MEM spatial eigenvectors...\n")

xy_mat   <- as.matrix(df_model[, c("x", "y")])
dist_mat <- dist(xy_mat)

# Build neighbour list and spatial weights
nb       <- spdep::dnearneigh(xy_mat, d1 = 0, d2 = 40)  # neighbours within 50m — tune as needed
listw    <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

# Compute MEMs — positive MEMs capture broad spatial patterns,
# negative MEMs capture fine-scale patterns
mem      <- adespatial::mem(listw)
mem_df   <- as.data.frame(mem)
names(mem_df) <- paste0("MEM", seq_len(ncol(mem_df)))

# Keep only significant MEMs (positive Moran's I)
mem_scores <- adespatial::moran.randtest(mem, listw, nrepet = 99)
sig_mems   <- which(mem_scores$pvalue < 0.05)
cat("Significant MEMs:", length(sig_mems), "\n")

mem_df_sig <- mem_df[, sig_mems, drop = FALSE]
mem_cols   <- names(mem_df_sig)

# Add MEMs to model dataframe
df_model <- bind_cols(df_model, mem_df_sig)

# --- Fit RF: species only vs species + MEM ---
results <- list()

for (response in spectral_cols) {
  
  cat("\n=== Response:", response, "===\n")
  
  # Species only
  rf_species <- ranger(
    formula    = as.formula(paste(response, "~", paste(species_cols, collapse = "+"))),
    data       = df_model,
    num.trees  = 1000,
    importance = "permutation",
    seed       = 42
  )
  cat("Species only  R²:", round(rf_species$r.squared, 3), "\n")
  
  # Species + spatial MEMs
  all_preds  <- c(species_cols, mem_cols)
  rf_spatial <- ranger(
    formula    = as.formula(paste(response, "~", paste(all_preds, collapse = "+"))),
    data       = df_model,
    num.trees  = 1000,
    importance = "permutation",
    seed       = 42
  )
  cat("Species + MEM R²:", round(rf_spatial$r.squared, 3), "\n")
  
  # Check residual spatial autocorrelation
  resid_species <- df_model[[response]] - rf_species$predictions
  resid_spatial <- df_model[[response]] - rf_spatial$predictions
  
  moran_before <- spdep::moran.test(resid_species, listw, zero.policy = TRUE)
  moran_after  <- spdep::moran.test(resid_spatial, listw, zero.policy = TRUE)
  
  cat("Moran's I before (species only):", round(moran_before$estimate[1], 4),
      " p =", round(moran_before$p.value, 4), "\n")
  cat("Moran's I after  (species + MEM):", round(moran_after$estimate[1], 4),
      " p =", round(moran_after$p.value, 4), "\n")
  
  # Variable importance — species only (exclude MEMs for interpretability)
  imp_df <- enframe(rf_spatial$variable.importance, name = "variable", value = "importance") %>%
    filter(variable %in% species_cols) %>%
    arrange(desc(importance)) %>%
    mutate(response = response)
  
  results[[response]] <- list(
    rf_species = rf_species,
    rf_spatial = rf_spatial,
    importance = imp_df,
    resid_before = resid_species,
    resid_after  = resid_spatial
  )
}

# --- R² comparison table ---
perf_df <- data.frame(
  response      = spectral_cols,
  r2_species    = sapply(spectral_cols, function(r) round(results[[r]]$rf_species$r.squared, 3)),
  r2_spatial    = sapply(spectral_cols, function(r) round(results[[r]]$rf_spatial$r.squared, 3))
)
cat("\n--- R² Comparison ---\n")
print(perf_df)

# --- Top 15 species importance plot (from spatial model) ---
imp_all <- bind_rows(lapply(results, `[[`, "importance")) %>%
  mutate(response = recode(response,
                           "global_PCA_20m_masked_5nm"     = "Global PCA",
                           "SA_entropy_20m_masked_5nm"     = "SA Entropy",
                           "SumBandEntropy_20m_masked_5nm" = "Sum Band Entropy"
  ))

top15 <- imp_all %>%
  group_by(response) %>%
  slice_max(importance, n = 15) %>%
  ungroup()

ggplot(top15, aes(x = reorder(variable, importance), y = importance, fill = response)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~response, scales = "free") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text  = element_text(face = "bold")
  ) +
  labs(
    x     = "Species",
    y     = "Permutation Importance",
    title = "Spatial RF (MEM): Species Composition vs Spectral Diversity"
  )

# --- Map residuals before and after spatial correction ---
for (response in spectral_cols) {
  plot_df <- centroids %>%
    mutate(
      resid_before = results[[response]]$resid_before,
      resid_after  = results[[response]]$resid_after
    ) %>%
    tidyr::pivot_longer(cols = c(resid_before, resid_after),
                        names_to = "model", values_to = "residual") %>%
    mutate(model = recode(model,
                          "resid_before" = "Species Only",
                          "resid_after"  = "Species + MEM"
    ))
  
  p <- ggplot(plot_df) +
    geom_sf(aes(color = residual), size = 1.2) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "Residual") +
    facet_wrap(~model) +
    theme_minimal() +
    ggtitle(paste("Residuals —", response))
  
  print(p)
}

##############################################################################

library(pdp)

for (response in spectral_cols) {
  
  # Use spatial model and pull top species (excluding MEMs)
  top_sp <- results[[response]]$importance %>%
    filter(variable %in% species_cols) %>%
    slice_max(importance, n = 1) %>%
    pull(variable)
  
  # Build training df with all predictors used in spatial model
  model_df <- df_model %>%
    dplyr::select(all_of(response), all_of(species_cols), all_of(mem_cols)) %>%
    rename(y = all_of(response))
  
  pd <- partial(
    results[[response]]$rf_spatial,
    pred.var = top_sp,
    train    = model_df
  )
  
  p <- autoplot(pd) +
    geom_hline(yintercept = mean(model_df$y, na.rm = TRUE),
               linetype = "dashed", color = "gray50") +
    labs(
      x     = paste(top_sp, "(stem count)"),
      y     = "Predicted spectral diversity (yhat)",
      title = paste(response, "~", top_sp)
    ) +
    theme_minimal()
  
  print(p)
}











