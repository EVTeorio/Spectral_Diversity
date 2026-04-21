# =============================================================================
# TOPOGRAPHY × SPECTRAL HETEROGENEITY → PHYLOGENETIC DIVERSITY
# Core Question: How does topography influence the relationship between
# spectral heterogeneity and species/phylogenetic diversity?
# Sub-question: How much ADDITIONAL variation does topography explain
# beyond spectral entropy alone?
# =============================================================================

library(dplyr)
library(ggplot2)
library(sf)
library(janitor)
library(spdep)
library(ggpubr)
library(vegan)       # for varpart()
library(car)         # for vif()
library(interactions) # for interaction plots (install if needed)

# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

df <- read.csv("Indices_SHPs/20m_spectral_sp.csv") %>%
  clean_names() %>%
  filter(complete.cases(.))

df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)

# --- Variable sets -----------------------------------------------------------
#Options
spectral_vars <- c("global_pca_20m_masked_5nm",
                   "sa_entropy_20m_masked_5nm",
                   "sa_entropy_smooth",
                   "sa_entropy_smooth_masked",
                   "sa_entropy_smooth_masked_711",
                   "sum_band_entropy_20m_masked_5nm"
)
entropy_var <- "sa_entropy_20m_masked_5nm"   # primary spectral predictor

topo_vars <- c(
  "dmnnt_v",    # topographic variable
  "avg_lvt",    # average elevation
  "elvtn_r"     # elevation range (topographic heterogeneity)
)

pd_vars <- c("afaith_pd", "faith_pd", "rao_pd", "shnnn_d")
response_var <- "rao_pd"

# --- Scale all predictors ----------------------------------------------------

all_preds <- c(entropy_var, topo_vars)

df_sf <- df_sf %>%
  mutate(across(all_of(all_preds),
                ~ scale(.) %>% as.numeric(),
                .names = "{.col}_sc"))

# Convenience scaled names
entropy_sc   <- paste0(entropy_var, "_sc")
topo_sc      <- paste0(topo_vars, "_sc")

# =============================================================================
# 2. NESTED MODEL COMPARISON
# Core: Does topography explain variation BEYOND spectral entropy?
# =============================================================================

# Model A: Spectral entropy only (baseline)
m_entropy <- lm(
  as.formula(paste(response_var, "~", entropy_sc)),
  data = df_sf
)

# Model B: Topography only
m_topo <- lm(
  as.formula(paste(response_var, "~", paste(topo_sc, collapse = " + "))),
  data = df_sf
)

# Model C: Full additive model (entropy + topo)
m_full <- lm(
  as.formula(paste(
    response_var, "~",
    paste(c(entropy_sc, topo_sc), collapse = " + ")
  )),
  data = df_sf
)

# --- Model fit summary -------------------------------------------------------

r2  <- c(
  Entropy_only = summary(m_entropy)$r.squared,
  Topo_only    = summary(m_topo)$r.squared,
  Full         = summary(m_full)$r.squared
)

adj_r2 <- c(
  Entropy_only = summary(m_entropy)$adj.r.squared,
  Topo_only    = summary(m_topo)$adj.r.squared,
  Full         = summary(m_full)$adj.r.squared
)

delta_r2_over_entropy <- r2["Full"] - r2["Entropy_only"]
delta_r2_over_topo    <- r2["Full"] - r2["Topo_only"]

cat("\n", strrep("=", 65), "\n")
cat("NESTED MODEL R² COMPARISON\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Entropy-only  R² = %.4f  |  Adj-R² = %.4f\n", r2[1], adj_r2[1]))
cat(sprintf("  Topo-only     R² = %.4f  |  Adj-R² = %.4f\n", r2[2], adj_r2[2]))
cat(sprintf("  Full model    R² = %.4f  |  Adj-R² = %.4f\n", r2[3], adj_r2[3]))
cat(sprintf("\n  ΔR² (Entropy → Full):  +%.4f  [topography's unique gain]\n",
            delta_r2_over_entropy))
cat(sprintf("  ΔR² (Topo → Full):     +%.4f  [entropy's unique gain]\n\n",
            delta_r2_over_topo))

# AIC comparison
aic_table <- AIC(m_entropy, m_topo, m_full)
cat("AIC Comparison:\n")
print(aic_table)

# Nested F-test: does adding topo significantly improve over entropy alone?
cat("\nNested F-test: Entropy-only vs Full model\n")
anova_result <- anova(m_entropy, m_full)
print(anova_result)

cat("\nNested F-test: Topo-only vs Full model\n")
print(anova(m_topo, m_full))

# VIF for full model (check multicollinearity)
cat("\nVariance Inflation Factors (Full model):\n")
print(vif(m_full))

# =============================================================================
# 3. VARIANCE PARTITIONING
# Cleanly separates: [entropy unique] [topo unique] [shared] [unexplained]
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("VARIANCE PARTITIONING (vegan::varpart)\n")
cat(strrep("=", 65), "\n")

# varpart uses matrix form — pull numeric data
Y  <- df_sf[[response_var]]
X1 <- as.matrix(df_sf[, entropy_sc])           # entropy fraction
X2 <- as.matrix(df_sf[, topo_sc])              # topo fraction

vp <- varpart(Y, X1, X2)
print(vp)

# Fractions:
# [a] = entropy unique
# [b] = shared (entropy ∩ topo)
# [c] = topo unique
# [d] = unexplained

cat("\nInterpretation:\n")
cat(sprintf(
  "  Entropy unique  [a]:  %.4f (%.1f%% of total explained)\n",
  vp$part$fract$Adj.R.square[1],
  vp$part$fract$Adj.R.square[1] / r2["Full"] * 100
))
cat(sprintf(
  "  Topo unique     [c]:  %.4f (%.1f%% of total explained)\n",
  vp$part$fract$Adj.R.square[3],
  vp$part$fract$Adj.R.square[3] / r2["Full"] * 100
))
cat(sprintf(
  "  Shared [b]:           %.4f (%.1f%% of total explained)\n",
  vp$part$fract$Adj.R.square[2],
  vp$part$fract$Adj.R.square[2] / r2["Full"] * 100
))

# Test significance of each fraction
cat("\nSignificance of entropy fraction (controlling for topo):\n")
print(anova(rda(Y, X1, X2)))

cat("\nSignificance of topo fraction (controlling for entropy):\n")
print(anova(rda(Y, X2, X1)))

# =============================================================================
# 4. INTERACTION MODELS
# Does topography MODERATE the entropy–PD relationship?
# i.e., does the slope of entropy→PD change across topo gradients?
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("INTERACTION MODELS: Does topo moderate entropy → PD?\n")
cat(strrep("=", 65), "\n")

# Test each topographic variable as a moderator separately to avoid collinearity
interaction_results <- lapply(topo_sc, function(tv) {
  form <- as.formula(paste(
    response_var, "~",
    entropy_sc, "*", tv
  ))
  m <- lm(form, data = df_sf)
  list(
    topo_var   = tv,
    model      = m,
    r2         = summary(m)$r.squared,
    adj_r2     = summary(m)$adj.r.squared,
    coefs      = summary(m)$coefficients,
    aic        = AIC(m)
  )
})

cat("\nInteraction model summaries (entropy × each topo var):\n\n")
for (res in interaction_results) {
  cat(strrep("-", 50), "\n")
  cat("Entropy ×", res$topo_var, "\n")
  cat(sprintf("  R² = %.4f  |  Adj-R² = %.4f  |  AIC = %.2f\n",
              res$r2, res$adj_r2, res$aic))
  int_row <- res$coefs[grep(":", rownames(res$coefs)), , drop = FALSE]
  if (nrow(int_row) > 0) {
    cat("  Interaction term(s):\n")
    print(round(int_row, 4))
  }
  cat("\n")
}

# Best-fit interaction model (lowest AIC)
best_int <- interaction_results[[which.min(sapply(interaction_results, `[[`, "aic"))]]
cat("Best interaction model (lowest AIC): entropy ×", best_int$topo_var, "\n")

# Full model vs best interaction model
cat("\nFull additive vs best interaction model:\n")
print(AIC(m_full, best_int$model))

# =============================================================================
# 5. SPATIAL AUTOCORRELATION CHECK
# Validate that residuals are independent before trusting p-values
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("SPATIAL AUTOCORRELATION: Moran's I on full model residuals\n")
cat(strrep("=", 65), "\n")

coords <- st_coordinates(st_centroid(df_sf))
nb     <- knn2nb(knearneigh(coords, k = 5))
lw     <- nb2listw(nb, style = "W")

moran_result <- moran.test(residuals(m_full), lw)
print(moran_result)

if (moran_result$p.value < 0.05) {
  cat("\n⚠  Significant spatial autocorrelation detected in residuals.\n")
  cat("   Consider spatial error or lag models (spatialreg package).\n")
  cat("   Current OLS R² estimates may be overconfident.\n\n")
  
  # Optional: simple spatial error model using spatialreg
  # library(spatialreg)
  # m_spatial <- errorsarlm(
  #   as.formula(paste(response_var, "~", paste(c(entropy_sc, topo_sc), collapse=" + "))),
  #   data = df_sf, listw = lw
  # )
  # summary(m_spatial)
} else {
  cat("\n✓  No significant spatial autocorrelation. OLS estimates are reliable.\n\n")
}

# =============================================================================
# 6. COEFFICIENT TABLE SUMMARY (publication-ready)
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("FULL MODEL COEFFICIENTS\n")
cat(strrep("=", 65), "\n")

full_coef <- as.data.frame(summary(m_full)$coefficients) %>%
  tibble::rownames_to_column("Term") %>%
  rename(
    Estimate   = Estimate,
    SE         = `Std. Error`,
    t_value    = `t value`,
    p_value    = `Pr(>|t|)`
  ) %>%
  mutate(
    Sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  )

print(full_coef, digits = 4)

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

# --- 7A. Variance partition pie (manual) -------------------------------------

vp_fractions <- data.frame(
  Component = c("Entropy unique", "Shared", "Topo unique", "Unexplained"),
  Proportion = c(
    max(0, vp$part$fract$Adj.R.square[1]),
    max(0, vp$part$fract$Adj.R.square[2]),
    max(0, vp$part$fract$Adj.R.square[3]),
    max(0, 1 - r2["Full"])
  )
)

p_vp <- ggplot(vp_fractions, aes(x = "", y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.8) +
  coord_polar("y") +
  scale_fill_manual(values = c(
    "Entropy unique" = "#2196F3",
    "Shared"         = "#9C27B0",
    "Topo unique"    = "#4CAF50",
    "Unexplained"    = "#E0E0E0"
  )) +
  geom_text(aes(label = ifelse(Proportion > 0.02,
                               paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  theme_void(base_size = 13) +
  labs(title = "Variance Partitioning: Rao PD",
       subtitle = "Spectral Entropy vs. Topography",
       fill = NULL) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

# --- 7B. Entropy → PD scatter, colored by elevation range -------------------

p_scatter <- ggplot(df_sf,
                    aes(x = .data[[entropy_sc]],
                        y = .data[[response_var]],
                        color = .data[[topo_sc[3]]])) +  # elvtn_r_sc
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE,
              linewidth = 1.2, linetype = "dashed") +
  scale_color_gradient2(low = "#1565C0", mid = "#FFF176",
                        high = "#B71C1C", midpoint = 0,
                        name = "Elev. Range\n(scaled)") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Spectral Entropy vs Rao PD",
    subtitle = "Points colored by topographic elevation range",
    x = "Spectral Entropy (scaled)",
    y = "Rao Phylogenetic Diversity"
  ) +
  theme(plot.title = element_text(face = "bold"))

# --- 7C. Model R² bar chart --------------------------------------------------

r2_df <- data.frame(
  Model = factor(c("Entropy only", "Topo only", "Full (additive)"),
                 levels = c("Entropy only", "Topo only", "Full (additive)")),
  R2    = as.numeric(r2)
)

p_r2 <- ggplot(r2_df, aes(x = Model, y = R2, fill = Model)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("R² = %.3f", R2)),
            vjust = -0.4, fontface = "bold", size = 4.2) +
  scale_fill_manual(values = c("#2196F3", "#4CAF50", "#9C27B0")) +
  scale_y_continuous(limits = c(0, max(r2) * 1.18), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Model R² Comparison",
    subtitle = paste0("ΔR² (Entropy→Full) = +", round(delta_r2_over_entropy, 3)),
    y = "R²", x = NULL
  ) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11))

# --- 7D. Partial residual plots (topo variables after entropy) ---------------

# Partial residuals of response after removing entropy
resid_after_entropy <- residuals(lm(
  as.formula(paste(response_var, "~", entropy_sc)), data = df_sf
))

df_partial <- df_sf %>%
  mutate(partial_resid = resid_after_entropy)

partial_plots <- lapply(seq_along(topo_sc), function(i) {
  tv  <- topo_sc[i]
  lbl <- c("Dominant Veg (dmnnt_v)",
           "Avg Live Veg (avg_lvt)",
           "Elevation Range (elvtn_r)")[i]
  
  ggplot(df_partial, aes(x = .data[[tv]], y = partial_resid)) +
    geom_point(alpha = 0.4, color = "#546E7A", size = 1.5) +
    geom_smooth(method = "lm", color = "#E53935", se = TRUE,
                linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    theme_minimal(base_size = 11) +
    labs(
      title = lbl,
      x = paste(lbl, "(scaled)"),
      y = "Residual Rao PD\n(after entropy)"
    ) +
    theme(plot.title = element_text(face = "bold", size = 10))
})

p_partial <- ggarrange(plotlist = partial_plots, ncol = 3, nrow = 1)
p_partial <- annotate_figure(p_partial,
                             top = text_grob(
                               "Partial Effects of Topography on Rao PD (after accounting for spectral entropy)",
                               face = "bold", size = 13
                             )
)

# --- Assemble & save all plots -----------------------------------------------

ggsave("vp_pie.png",          plot = p_vp,      width = 6,  height = 5,  dpi = 300)
ggsave("entropy_pd_scatter.png", plot = p_scatter, width = 7,  height = 5,  dpi = 300)
ggsave("model_r2_bars.png",   plot = p_r2,      width = 6,  height = 4.5, dpi = 300)
ggsave("partial_topo.png",    plot = p_partial,  width = 12, height = 4,  dpi = 300)

cat("\n✓  Plots saved: vp_pie.png | entropy_pd_scatter.png | model_r2_bars.png | partial_topo.png\n")

# =============================================================================
# 8. FINAL NARRATIVE SUMMARY
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("ECOLOGICAL INTERPRETATION SUMMARY\n")
cat(strrep("=", 65), "\n")

cat(sprintf("
STEP 1 — BASELINE (Entropy → Rao PD):
  Spectral entropy alone explains %.1f%% of variation in Rao PD.
  This is your baseline spectral heterogeneity signal.
 
STEP 2 — TOPOGRAPHIC ADDITION:
  Adding topographic variables increases R² by %.1f percentage points
  (from %.3f to %.3f).
 
STEP 3 — UNIQUE CONTRIBUTIONS (Variance Partitioning):
  - Spectral entropy unique:  see [a] fraction above
  - Topography unique:        see [c] fraction above
  - Shared/collinear portion: see [b] fraction above
  (If shared is large, topo and entropy are capturing overlapping gradients)
 
STEP 4 — MODERATION TEST:
  Interaction models test whether the slope of entropy → PD
  changes depending on topographic context (e.g., steeper entropy
  effect in high-relief terrain). See interaction results above.
 
STEP 5 — SPATIAL VALIDITY:
  Check Moran's I result above. If autocorrelation is present,
  treat p-values and R² with caution until spatial models are fit.
",
            r2["Entropy_only"] * 100,
            delta_r2_over_entropy * 100,
            r2["Entropy_only"],
            r2["Full"]
))
