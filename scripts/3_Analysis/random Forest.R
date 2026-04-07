
library(dplyr)
library(ranger)
library(ggplot2)
library(tidyr)
library(tibble)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
df <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  filter(complete.cases(.))

head(df)# --- Setup ---
spectral_cols <- c("global_PCA_20m_masked_5nm", "SA_entropy_20m_masked_5nm", 
                   "SumBandEntropy_20m_masked_5nm")

species_start <- which(names(df) == "CAGL8")
species_cols  <- names(df)[species_start:ncol(df)]

# Drop rows with NA in any spectral or species column
df_clean <- df %>%
  dplyr::select(Name, all_of(spectral_cols), all_of(species_cols)) %>%
  drop_na()

# --- Ranger RF for each spectral diversity response ---
results <- list()

for (response in spectral_cols) {
  
  cat("\n=== Response:", response, "===\n")
  
  model_df <- df_clean %>%
    dplyr::select(all_of(response), all_of(species_cols)) %>%
    rename(y = all_of(response))
  
  # Fit random forest
  rf <- ranger(
    formula       = y ~ .,
    data          = model_df,
    num.trees     = 1000,
    importance    = "permutation",
    respect.unordered.factors = "partition",
    seed          = 42
  )
  
  cat("R-squared (OOB):", round(rf$r.squared, 3), "\n")
  cat("RMSE (OOB):     ", round(sqrt(rf$prediction.error), 4), "\n")
  
  # Variable importance
  imp_df <- enframe(rf$variable.importance, name = "species", value = "importance") %>%
    arrange(desc(importance)) %>%
    mutate(response = response)
  
  results[[response]] <- list(rf = rf, importance = imp_df)
}

# --- Plot variable importance: top 15 species per response ---
imp_all <- bind_rows(lapply(results, `[[`, "importance"))

# Clean up response names for plot labels
imp_all <- imp_all %>%
  mutate(response = recode(response,
                           "global_PCA_20m_masked_5nm"    = "Global PCA",
                           "SA_entropy_20m_masked_5nm"    = "SA Entropy",
                           "SumBandEntropy_20m_masked_5nm"= "Sum Band Entropy"
  ))

# Top 15 per response
top15 <- imp_all %>%
  group_by(response) %>%
  slice_max(importance, n = 15) %>%
  ungroup()

ggplot(top15, aes(x = reorder(species, importance), y = importance, fill = response)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~response, scales = "free") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 9),
    strip.text   = element_text(face = "bold")
  ) +
  labs(
    x     = "Species",
    y     = "Permutation Importance",
    title = "Random Forest: Species Composition vs Spectral Diversity"
  )

# --- Summary table of model performance ---
perf_df <- data.frame(
  response  = spectral_cols,
  r_squared = sapply(spectral_cols, function(r) round(results[[r]]$rf$r.squared, 3)),
  rmse      = sapply(spectral_cols, function(r) round(sqrt(results[[r]]$rf$prediction.error), 4))
)

cat("\n--- Model Performance Summary ---\n")
print(perf_df)

# --- Partial dependence for top species per response (optional deep dive) ---
# Uncomment to run — shows how each top species drives the response
library(pdp)
for (response in spectral_cols) {
  top_sp <- results[[response]]$importance$species[4]
  model_df <- df_clean %>%
    dplyr::select(all_of(response), all_of(species_cols)) %>%
    rename(y = all_of(response))
  pd <- partial(results[[response]]$rf, pred.var = top_sp, train = model_df)
  print(autoplot(pd) + ggtitle(paste(response, "~", top_sp)))
}
