
# =============================================================================
# Spectral + Diversity Analysis Report
# =============================================================================

library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(nlme)
library(janitor)
library(stringr)
library(readr)
library(gridExtra)
library(grid)
library(ggpubr)


setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
df <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  clean_names() %>%
  filter(complete.cases(.))

# ---------------------------
# 2) Convert geometry to sf
# ---------------------------
df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)

# Compute point-on-surface coordinates
pt <- st_point_on_surface(df_sf$geometry)
xy <- st_coordinates(pt)
df_sf <- df_sf %>%
  mutate(lon_mean = xy[, 1],
         lat_mean = xy[, 2])

# ---------------------------
# 3) Define predictor & response variables
# ---------------------------
spectral_vars <- c("global_pca_20m_masked_5nm",
                   "sa_entropy_20m_masked_5nm",
                   "sum_band_entropy_20m_masked_5nm")

pd_vars <- c("afaith_pd", "faith_pd", "rao_pd", "rchnss_", "shnnn_d", "smpsn_d")

# For demonstration, pick one PD metric as response (e.g., faith_pd)
response_var <- "afaith_pd"
predictor_var <- "global_pca_20m_masked_5nm"

# ---------------------------
# 4) Scale predictors
# ---------------------------
df_sf <- df_sf %>%
  mutate(across(all_of(spectral_vars), ~scale(.) %>% as.numeric(), .names = "{.col}_sc"))
# df_sf <- df_sf %>%
#   mutate(across(all_of(spectral_vars), as.numeric))

# Save means and SDs
means <- sapply(df_sf[spectral_vars], mean, na.rm = TRUE)
sds   <- sapply(df_sf[spectral_vars], sd, na.rm = TRUE)

# ---------------------------
# 5) Exploratory plots
# ---------------------------
for(var in spectral_vars) {
  p <- ggplot(df_sf, aes_string(x = var)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black") +
    labs(title = paste("Distribution of", var), x = var, y = "Count") +
    theme_minimal()
  print(p)
}

ggplot(df_sf, aes_string(x = predictor_var, y = response_var)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  labs(title = paste(response_var, "vs", predictor_var)) +
  theme_minimal()

# ---------------------------
# 6) Correlations
# ---------------------------
cor(df_sf[[response_var]], df_sf[[predictor_var]], method = "pearson")
cor(df_sf[[response_var]], df_sf[[predictor_var]], method = "spearman")

# ---------------------------
# 7) Ordinary Least Squares (OLS)
# ---------------------------
formula_ols <- as.formula(paste(response_var, "~", paste0(predictor_var, "_sc")))
ols_model <- lm(formula_ols, data = df_sf)
summary(ols_model)

# Residuals and fitted values
df_sf <- df_sf %>%
  mutate(fitted_ols = fitted(ols_model),
         resid_ols  = resid(ols_model))

# ---------------------------
# 8) Spatial coordinates in meters (UTM projection)
# ---------------------------
pts_ll <- st_as_sf(df_sf, coords = c("lon_mean", "lat_mean"), crs = 4326)
pts_utm <- st_transform(pts_ll, 26916)
centroids <- st_centroid(df_sf)
xy_m <- st_coordinates(centroids)

df_sf <- df_sf %>%
  mutate(x_m = xy_m[,1],
         y_m = xy_m[,2])

# ---------------------------
# 9) Nearest-neighbor weights & Moran's I
# ---------------------------
knn <- knearneigh(xy_m, k = 8)
nb  <- knn2nb(knn)
lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran.test(df_sf[[response_var]], lw, zero.policy = TRUE)
moran.test(df_sf$resid_ols, lw, zero.policy = TRUE)

# ---------------------------
# 10) Spatial GLS
# ---------------------------
gls_model <- gls(
  formula_ols,
  data = df_sf,
  correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
  method = "REML"
)

summary(gls_model)

df_sf <- df_sf %>%
  mutate(fitted_gls = predict(gls_model),
         resid_gls  = residuals(gls_model, type = "normalized"))

moran.test(df_sf$resid_gls, lw, zero.policy = TRUE)

# ---------------------------
# 11) Model comparison (ML)
# ---------------------------
gls_null <- gls(as.formula(paste(response_var, "~1")),
                data = df_sf,
                correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
                method = "ML")

gls_predictor <- gls(formula_ols,
                     data = df_sf,
                     correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
                     method = "ML")

AIC(gls_null, gls_predictor)

# ---------------------------
# 12) Likelihood-ratio pseudo-R2
# ---------------------------
n <- nrow(df_sf)
ll0 <- as.numeric(logLik(gls_null))
ll1 <- as.numeric(logLik(gls_predictor))

R2_LR <- 1 - exp(-(2*(ll1 - ll0))/n)
R2_LR

# ---------------------------
# 13) Predictions
# ---------------------------

newdat <- setNames(
  data.frame(c(0, 1, -1)),
  paste0(predictor_var, "_sc")
)

newdat$pred_ols <- predict(ols_model, newdata = newdat)
newdat$pred_gls <- predict(gls_model, newdata = newdat)
newdat
#################################################################################
#            Report Out                 ####################################
##############################################################################
sig_label <- function(p) {
  if (p < 0.001) return("highly significant (p < 0.001)")
  if (p < 0.01)  return("significant (p < 0.01)")
  if (p < 0.05)  return("significant (p < 0.05)")
  return("not significant (p >= 0.05)")
}

direction <- ifelse(r_pearson > 0, "positively", "negatively")

plain_language <- paste0(
  "=============================================================
PLAIN LANGUAGE SUMMARY
=============================================================

VARIABLES ANALYZED
  Response  : ", response_var, "
  Predictor : ", predictor_var, "
  Sample size: ", n, " quads

-------------------------------------------------------------
CORRELATION
-------------------------------------------------------------
  The response variable (", response_var, ") and predictor
  (", predictor_var, ") are ", direction, " correlated.
  Pearson r  = ", round(r_pearson, 3), "  (p = ",
  round(cor_test$p.value, 4), " — ", sig_label(cor_test$p.value), ")
  Spearman r = ", round(r_spearman, 3), "
  Interpretation: A Pearson r of ", round(r_pearson, 3), " means that as
  ", predictor_var, " increases, ", response_var, " tends to ",
  ifelse(r_pearson > 0, "increase", "decrease"), ".

-------------------------------------------------------------
OLS REGRESSION (ignoring spatial structure)
-------------------------------------------------------------
  R² = ", round(ols_r2, 3), "
  This means ", round(ols_r2 * 100, 1), "% of the variation in ",
  response_var, " is explained by ", predictor_var, "
  before accounting for spatial autocorrelation.

  Model F-statistic = ", round(ols_fstat[1], 3),
  "  (p = ", round(ols_pval, 4), " — ", sig_label(ols_pval), ")
  Coefficient for ", predictor_var, " (scaled):
    Estimate = ", round(ols_coef[2, 1], 4), "
    Std. Error = ", round(ols_coef[2, 2], 4), "
    p-value = ", round(ols_coef[2, 4], 4), " — ",
  sig_label(ols_coef[2, 4]), "

-------------------------------------------------------------
SPATIAL AUTOCORRELATION CHECK
-------------------------------------------------------------
  Moran's I of ", response_var, " itself:
    I = ", round(moran_response$estimate[1], 4),
  "  (p = ", round(moran_response$p.value, 4), " — ",
  sig_label(moran_response$p.value), ")
  Interpretation: ", ifelse(moran_response$p.value < 0.05,
                            paste0("Significant spatial autocorrelation is present in ",
                                   response_var, ". Nearby quads tend to have similar values,
  which violates OLS independence assumptions and inflates
  Type I error. A spatial model is warranted."),
                            paste0("No significant spatial autocorrelation detected in ",
                                   response_var, ". OLS assumptions may be adequate.")), "

  Moran's I of OLS residuals:
    I = ", round(moran_resid_ols$estimate[1], 4),
  "  (p = ", round(moran_resid_ols$p.value, 4), " — ",
  sig_label(moran_resid_ols$p.value), ")
  Interpretation: ", ifelse(moran_resid_ols$p.value < 0.05,
                            "Spatial autocorrelation remains in OLS residuals, confirming
  that OLS results are unreliable for inference here.",
                            "OLS residuals show no significant spatial autocorrelation."), "

  Moran's I of GLS residuals (after spatial correction):
    I = ", round(moran_resid_gls$estimate[1], 4),
  "  (p = ", round(moran_resid_gls$p.value, 4), " — ",
  sig_label(moran_resid_gls$p.value), ")
  Interpretation: ", ifelse(moran_resid_gls$p.value < 0.05,
                            "Some spatial autocorrelation remains after GLS correction.
  Consider a more flexible spatial covariance structure.",
                            "GLS successfully removed spatial autocorrelation from residuals.
  GLS inference is reliable."), "

-------------------------------------------------------------
SPATIAL GLS (accounting for spatial autocorrelation)
-------------------------------------------------------------
  Likelihood-ratio pseudo R² = ", round(R2_LR, 3), "
  This means approximately ", round(R2_LR * 100, 1), "% of variance in
  ", response_var, " is explained by ", predictor_var, "
  after accounting for spatial structure.

  Likelihood-ratio test vs null spatial model:
    LRT statistic = ", round(lrt_stat, 3), "
    p-value = ", round(lrt_pval, 4), " — ", sig_label(lrt_pval), "

  GLS coefficient for ", predictor_var, " (scaled):
    Estimate  = ", round(gls_sum$tTable[2, 1], 4), "
    Std. Error = ", round(gls_sum$tTable[2, 2], 4), "
    p-value   = ", round(gls_sum$tTable[2, 4], 4), " — ",
  sig_label(gls_sum$tTable[2, 4]), "

  AIC comparison:
    Null spatial model AIC  = ", round(aic_table[1, 2], 2), "
    Predictor model AIC     = ", round(aic_table[2, 2], 2), "
    ΔAIC = ", round(aic_table[1, 2] - aic_table[2, 2], 2), "
  Interpretation: ", ifelse((aic_table[1, 2] - aic_table[2, 2]) > 2,
                            paste0("The predictor model is better supported (ΔAIC > 2), confirming
  that ", predictor_var, " explains meaningful variance in ",
                                   response_var, " even after spatial correction."),
                            paste0("The predictor model offers little improvement over the null
  spatial model (ΔAIC <= 2). The relationship may be largely
  driven by shared spatial structure rather than a direct
  ecological association.")), "

-------------------------------------------------------------
OVERALL INTERPRETATION
-------------------------------------------------------------
  Before spatial correction: OLS R² = ", round(ols_r2, 3), "
  After  spatial correction: pseudo R² = ", round(R2_LR, 3), "
  Change in variance explained: ",
  round((ols_r2 - R2_LR) * 100, 1), " percentage points

  ", ifelse(ols_r2 > R2_LR,
            paste0("The drop in R² after spatial correction (",
                   round((ols_r2 - R2_LR) * 100, 1),
                   " pp) suggests that part of the apparent relationship
  between ", predictor_var, " and ", response_var, "
  was due to shared spatial structure (both variables
  being spatially clustered) rather than a direct
  ecological effect."),
            "R² was maintained or increased after spatial correction,
  suggesting the relationship is robust to spatial structure."), "
=============================================================
"
)

cat(plain_language)

# ---------------------------
# 14) Compile PDF report
# ---------------------------
pdf_path <- paste0(
  "Spectral_Diversity/Documents/",
  response_var, "_vs_", predictor_var, "_report.pdf"
)

summary_grob <- textGrob(
  plain_language,
  x = 0.05, y = 0.95,
  hjust = 0, vjust = 1,
  gp = gpar(fontsize = 7, fontfamily = "mono")
)

pdf(pdf_path, width = 11, height = 8.5)

# Page 1 — plain language summary
grid.newpage()
grid.draw(summary_grob)

# Page 2 — histograms
grid.newpage()
do.call(grid.arrange, c(hist_plots, ncol = 3,
                        top = "Distributions of Spectral Diversity Variables"))

# Page 3 — scatter + residual plot
grid.newpage()
grid.arrange(scatter_plot, resid_plot_ols, ncol = 2,
             top = "OLS: Scatter Plot and Residual Diagnostics")

# Page 4 — Moran plot + OLS residual map
grid.newpage()
grid.arrange(moran_plot, resid_map_ols, ncol = 2,
             top = "Spatial Autocorrelation Diagnostics")

# Page 5 — GLS residual map
grid.newpage()
grid.arrange(resid_map_gls, ncol = 1,
             top = "GLS Residuals After Spatial Correction")

dev.off()

cat("\nPDF report written to:", pdf_path, "\n")