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
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
df <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  clean_names() %>%
  dplyr::filter(complete.cases(.))

# ---------------------------
# 2) Convert geometry to sf
# ---------------------------
df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)

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
                   "sa_entropy_smooth",
                   "sa_entropy_smooth_masked",
                   "sa_entropy_smooth_masked_711",
                   "sum_band_entropy_20m_masked_5nm"
                   )

pd_vars <- c("afaith_pd", "faith_pd", "rao_pd", "rchnss_", "shnnn_d", "smpsn_d")

response_var  <- "afaith_pd"
predictor_var <- "sa_entropy_smooth_masked_711"

# ---------------------------
# 4) Scale predictors
# ---------------------------
df_sf <- df_sf %>%
  mutate(across(all_of(spectral_vars), ~scale(.) %>% as.numeric(),
                .names = "{.col}_sc"))

means <- sapply(df_sf[spectral_vars], mean, na.rm = TRUE)
sds   <- sapply(df_sf[spectral_vars], sd, na.rm = TRUE)

# ---------------------------
# 5) Exploratory plots
# ---------------------------
hist_plots <- list()
for (var in spectral_vars) {
  p <- ggplot(df_sf, aes_string(x = var)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black") +
    labs(title = paste("Distribution of", var), x = var, y = "Count") +
    theme_minimal()
  hist_plots[[var]] <- p
}

scatter_plot <- ggplot(df_sf, aes_string(x = predictor_var, y = response_var)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(
    title    = paste(response_var, "vs", predictor_var),
    subtitle = "Red line = OLS regression fit",
    x        = predictor_var,
    y        = response_var
  ) +
  theme_minimal()

# ---------------------------
# 6) Correlations
# ---------------------------
r_pearson  <- cor(df_sf[[response_var]], df_sf[[predictor_var]], method = "pearson")
r_spearman <- cor(df_sf[[response_var]], df_sf[[predictor_var]], method = "spearman")

cor_test   <- cor.test(df_sf[[response_var]], df_sf[[predictor_var]], method = "pearson")

# ---------------------------
# 7) OLS
# ---------------------------
formula_ols <- as.formula(paste(response_var, "~", paste0(predictor_var, "_sc")))
ols_model   <- lm(formula_ols, data = df_sf)
ols_summary <- summary(ols_model)

df_sf <- df_sf %>%
  mutate(fitted_ols = fitted(ols_model),
         resid_ols  = resid(ols_model))

ols_r2    <- ols_summary$r.squared
ols_coef  <- ols_summary$coefficients
ols_fstat <- ols_summary$fstatistic
ols_pval  <- pf(ols_fstat[1], ols_fstat[2], ols_fstat[3], lower.tail = FALSE)

# Residual plot
resid_plot_ols <- ggplot(df_sf, aes(x = fitted_ols, y = resid_ols)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
  labs(title = "OLS Residuals vs Fitted",
       x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# ---------------------------
# 8) Spatial coordinates
# ---------------------------
pts_utm  <- st_transform(
  st_as_sf(df_sf, coords = c("lon_mean", "lat_mean"), crs = 4326), 26916)
centroids <- st_centroid(df_sf)
xy_m      <- st_coordinates(centroids)

df_sf <- df_sf %>%
  mutate(x_m = xy_m[, 1],
         y_m = xy_m[, 2])

# ---------------------------
# 9) Moran's I
# ---------------------------
knn <- knearneigh(xy_m, k = 8)
nb  <- knn2nb(knn)
lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_response <- moran.test(df_sf[[response_var]], lw, zero.policy = TRUE)
moran_resid_ols <- moran.test(df_sf$resid_ols, lw, zero.policy = TRUE)

# Moran plot
moran_plot <- ggplot(df_sf, aes(x = .data[[response_var]],
                                y = lag(df_sf[[response_var]]))) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(title = paste("Moran Plot —", response_var),
       x = response_var,
       y = paste("Spatially Lagged", response_var)) +
  theme_minimal()

# Residual spatial map
resid_map_ols <- ggplot(df_sf) +
  geom_sf(aes(color = resid_ols), size = 1.2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "Residual") +
  labs(title = "Spatial Distribution of OLS Residuals") +
  theme_minimal()

# ---------------------------
# 10) Spatial GLS
# ---------------------------
gls_model <- gls(
  formula_ols,
  data        = df_sf,
  correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
  method      = "REML"
)

gls_sum <- summary(gls_model)

df_sf <- df_sf %>%
  mutate(fitted_gls = predict(gls_model),
         resid_gls  = residuals(gls_model, type = "normalized"))

moran_resid_gls <- moran.test(df_sf$resid_gls, lw, zero.policy = TRUE)

gls_scatter <- ggplot(df_sf, aes(x = .data[[predictor_var]], y = .data[[response_var]])) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_line(aes(y = fitted_gls), color = "firebrick", linewidth = 1.1) +
  labs(
    title    = paste("GLS:", response_var, "vs", predictor_var),
    subtitle = "Red line = GLS fitted values (spatial autocorrelation corrected)",
    x        = predictor_var,
    y        = response_var
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )


print(gls_scatter)



# ---------------------------
# 11) Model comparison
# ---------------------------
gls_null <- gls(as.formula(paste(response_var, "~ 1")),
                data        = df_sf,
                correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
                method      = "ML")

gls_pred <- gls(formula_ols,
                data        = df_sf,
                correlation = corExp(form = ~ x_m + y_m, nugget = TRUE),
                method      = "ML")

aic_table <- AIC(gls_null, gls_pred)

# ---------------------------
# 12) Likelihood-ratio pseudo R2
# ---------------------------
n    <- nrow(df_sf)
ll0  <- as.numeric(logLik(gls_null))
ll1  <- as.numeric(logLik(gls_pred))
R2_LR <- 1 - exp(-(2 * (ll1 - ll0)) / n)

lrt_stat <- 2 * (ll1 - ll0)
lrt_pval <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# ---------------------------
# 13) Plain language summary
# ---------------------------
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
                            paste0("Significant spatial autocorrelation in ",
                                   response_var, ".Nearby quads tend to have similar values,
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
  Likelihood-ratio pseudo R² = ", round(R2_LR, 3),"
  This means ", round(R2_LR * 100, 1), "% of variance in
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
  that ", predictor_var, " explains additional variation in ",
                                   response_var, " beyond spatial structure."),
                            paste0("The predictor model offers little improvement over the null
  spatial model (ΔAIC <= 2). 
  The relationship may be largely
  driven by shared spatial structure.")), "

-------------------------------------------------------------
OVERALL INTERPRETATION
-------------------------------------------------------------
  Before spatial correction: OLS R² = ", round(ols_r2, 3), "
  After  spatial correction: pseudo R² = ", round(R2_LR, 3), "
  Change in variance explained: ",
  round((ols_r2 - R2_LR) * 100, 1), " percentage points

  ", ifelse(abs(ols_r2 - R2_LR) < 0.005,
            paste0("There is little to no change in explanatory power after spatial correction (ΔR² ≈ 0). 
  This suggests that the relationship between ",
                   predictor_var, " and ", response_var, "
  is robust to spatial autocorrelation and not primarily driven by shared spatial structure."),
            
            ifelse(ols_r2 > R2_LR,
                   paste0("The reduction in R² after spatial correction (",
                          round((ols_r2 - R2_LR) * 100, 1),
                          " pp) suggests that part of the apparent relationship
  between ", predictor_var, " and ", response_var, "
  was due to shared spatial structure (i.e., both variables being spatially clustered)."),
                   
                   paste0("The increase in explanatory power after spatial correction suggests that
  accounting for spatial structure clarifies the relationship between ",
                          predictor_var, " and ", response_var, "."))), "
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

print_text_pages <- function(text, fontsize = 12, fontfamily = "mono",
                             lines_per_page = 50) {
  all_lines <- unlist(strsplit(text, "\n"))
  pages     <- split(all_lines, ceiling(seq_along(all_lines) / lines_per_page))
  
  for (page_lines in pages) {
    grid.newpage()
    # Define a viewport with margins
    pushViewport(viewport(
      x      = 0.5,
      y      = 0.5,
      width  = 0.95,   # left/right margins
      height = 0.92,   # top/bottom margins
      just   = c("center", "center")
    ))
    grid.text(
      label  = paste(page_lines, collapse = "\n"),
      x      = 0,
      y      = 1,
      hjust  = 0,
      vjust  = 1,
      gp     = gpar(fontsize = fontsize, fontfamily = fontfamily)
    )
    popViewport()
  }
}

pdf(pdf_path, width = 8.5, height = 11)

# Page(s) 1+ — plain language summary
print_text_pages(plain_language, fontsize = 9, lines_per_page = 50)

# Plot pages
do.call(grid.arrange, c(hist_plots, ncol = 2,
                        top = "Distributions of Spectral Diversity Variables"))

grid.arrange(scatter_plot, resid_plot_ols, ncol = 1,
             top = "OLS: Scatter Plot and Residual Diagnostics")

grid.arrange(moran_plot, ncol = 1,
             top = "Spatial Autocorrelation Diagnostics")

dev.off()
cat("\nPDF report written to:", pdf_path, "\n")
beep()

