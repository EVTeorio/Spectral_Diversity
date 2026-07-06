# ============================================================
# SPATIAL ECOLOGICAL MODELING PIPELINE
# Spectral heterogeneity → phylogenetic diversity
# Includes: LM, GAMM (spatial), Random Forest
# ============================================================

library(dplyr)
library(ggplot2)
library(sf)
library(mgcv)
library(nlme)
library(janitor)
library(car)
library(AICcmodavg)
library(ranger)

# ---------------------------
# 1) DATA PREP
# ---------------------------

edge_quads <- c(
  "0","1","100","2","200","3","300","4","400","5","500","6","600",
  "7","700","8","800","9","900","10","1000","11","1100","12","1200","13",
  "14","1300","15","1400","16","1500","17","1600","18","1700","19","1800","20",
  "1900","21","22","1901","23","24","1903","124","1905","224","1906","324","1907",
  "424","1908","524","1909","624","1910","1911","724","1912","824","1913","924","1914",
  "1024","1915","1124","1916","1224","1324","1917","1424","1919","1524","1920","1624","1921",
  "1724","1922","1824","1923","1924","1904","1902","1918"
)

df_sf <- read.csv("Spectral_Diversity/Quad_Values/20m_spectral_sp.csv") %>%
  clean_names() %>%
  filter(complete.cases(.),
         !name %in% edge_quads)

df_sf$name <- as.character(df_sf$name)

# ---------------------------
# 2) VARIABLES
# ---------------------------

response_var <- "afaith_pd"

spectral_vars <- c(
  "sa_entropy_smooth_masked_711",
  "sa_entropy_smooth_masked",
  "sa_entropy_smooth",
  "global_pca_20m_masked_5nm",
  "sa_entropy_20m_masked_5nm",
  "sum_band_entropy_20m_masked_5nm"
)

topo_vars_cont <- c("avg_lvt", "elvtn_r")
topo_var_cat   <- "dmnnt_v"

# spatial coordinates from geometry
df_sf <- st_as_sf(df_sf, wkt = "geometry", crs = 4326)

cent <- st_centroid(df_sf)
xy <- st_coordinates(cent)

df_sf$x_m <- xy[,1]
df_sf$y_m <- xy[,2]

df_sf$dmnnt_v <- as.factor(df_sf$dmnnt_v)

# ---------------------------
# 3) SCALE VARIABLES
# ---------------------------

df_sf <- df_sf %>%
  mutate(
    across(all_of(c(spectral_vars, topo_vars_cont)),
           ~as.numeric(scale(.)),
           .names = "{.col}_sc")
  )

predictor_var <- "sa_entropy_20m_masked_5nm"

# ---------------------------
# 4) LINEAR MODELS
# ---------------------------

m_lm <- lm(
  as.formula(paste(response_var, "~", paste0(predictor_var, "_sc"))),
  data = df_sf
)

m_lm_topo <- lm(
  as.formula(paste(
    response_var, "~",
    paste0(predictor_var, "_sc"),
    "+ avg_lvt_sc + elvtn_r_sc + dmnnt_v"
  )),
  data = df_sf
)

m_inter <- lm(
  as.formula(paste(
    response_var, "~",
    paste0(predictor_var, "_sc"),
    "* dmnnt_v + avg_lvt_sc + elvtn_r_sc"
  )),
  data = df_sf
)

# ---------------------------
# 5) GAM (NON-SPATIAL)
# ---------------------------

m_gam <- gam(
  as.formula(paste(
    response_var, "~",
    paste0("s(", predictor_var, "_sc)"),
    "+ s(avg_lvt_sc) + s(elvtn_r_sc) + dmnnt_v"
  )),
  data = df_sf,
  method = "REML"
)

# ---------------------------
# 6) SPATIALLY CORRECTED GAMM
# ---------------------------

gamm_formula <- as.formula(
  paste(
    response_var, "~",
    paste0("s(", predictor_var, "_sc)"),
    "+ s(avg_lvt_sc) + s(elvtn_r_sc) + dmnnt_v"
  )
)

m_gamm <- gamm(
  formula = gamm_formula,
  data = df_sf,
  correlation = corExp(
    form = ~ x_m + y_m,
    nugget = TRUE
  ),
  method = "REML"
)

# ---------------------------
# 7) RANDOM FOREST (RANGER)
# ---------------------------

rf_formula <- as.formula(
  paste(
    response_var, "~",
    paste(
      c(
        paste0(spectral_vars, "_sc"),
        paste0(topo_vars_cont, "_sc"),
        "dmnnt_v"
      ),
      collapse = " + "
    )
  )
)

rf_data <- df_sf %>%
  st_drop_geometry()

rf_model <- ranger(
  formula = rf_formula,
  data = rf_data,
  importance = "impurity",
  num.trees = 500
)

df_sf$rf_pred <- predict(rf_model, data = rf_data)$predictions

# ---------------------------
# 8) MODEL COMPARISON
# ---------------------------

# GAM comparison (NO GAMM here)
aic_table <- AIC(
  m_lm,
  m_lm_topo,
  m_inter,
  m_gam
)

print(aic_table)

# spatial GAMM must be evaluated separately
cat("\n--- GAMM (spatial model) summary ---\n")
print(summary(m_gamm$lme))

anova(m_lm, m_lm_topo, m_inter)

# ---------------------------
# 9) GAMM OUTPUT
# ---------------------------

summary(m_gamm$gam)

# ---------------------------
# 10) RF IMPORTANCE
# ---------------------------

importance <- as.data.frame(rf_model$variable.importance)

# ---------------------------
# 11) VISUAL CHECK
# ---------------------------

df_sf$fitted_inter <- fitted(m_inter)

ggplot(df_sf, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = dmnnt_v), alpha = 0.6) +
  geom_point(aes(y = rf_pred), color = "darkgreen", alpha = 0.4) +
  geom_point(aes(y = fitted_inter), color = "red", alpha = 0.4) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral + Topography",
    subtitle = "Black = linear | Red = interaction | Green = RF predictions"
  )

# ---------------------------
# 12) RESIDUAL TEST (topo after spectral)
# ---------------------------

df_sf$resid_lm <- resid(m_lm)

lm_topo_resid <- lm(resid_lm ~ avg_lvt_sc + elvtn_r_sc + dmnnt_v, data = df_sf)
summary(lm_topo_resid)

# ---------------------------
# 13) SUMMARY INTERPRETATION
# ---------------------------

cat("
============================================================
MODEL SUMMARY
============================================================

1. Linear spectral R2:
", summary(m_lm)$r.squared, "

2. Full model R2:
", summary(m_lm_topo)$r.squared, "

3. GAM indicates nonlinear structure:
Check edf values in summary(m_gam)

4. GAMM accounts for spatial autocorrelation:
Look at correlation structure in residuals

5. Random Forest:
Captures nonlinear + interaction effects without assumptions

============================================================
")
