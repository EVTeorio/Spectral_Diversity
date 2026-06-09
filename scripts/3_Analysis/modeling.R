

library(dplyr)
library(ggplot2)
library(sf)
library(mgcv)
library(nlme)
library(janitor)
library(car)
library(AICcmodavg)

# ---------------------------
# 1) Data prep
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

df_sf <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  clean_names() %>%
  dplyr::filter(
    complete.cases(.),
    !name %in% edge_quads   # <-- comment this line to see in scatter plot.
  ) %>%
  st_as_sf(wkt = "geometry", crs = 4326)

# ensure numeric ID consistency
df_sf$name <- as.character(df_sf$name)

# ---------------------------
# 2) Define variables
# ---------------------------

response_vars <- c("afaith_pd", "faith_pd", "rao_pd", "shnnn_d")

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

response_var  <- "afaith_pd"
predictor_var <- "sa_entropy_20m_masked_5nm"

# ---------------------------
# 3) Scale predictors
# ---------------------------

df_sf <- df_sf %>%
  mutate(
    across(all_of(c(spectral_vars, topo_vars_cont)),
           ~as.numeric(scale(.)),
           .names = "{.col}_sc"),
    dmnnt_v = as.factor(dmnnt_v)
  )

# ---------------------------
# 4) MODEL SET
# ---------------------------

# (A) Linear model (spectral only)
m_lm <- lm(
  as.formula(paste(response_var, "~", paste0(predictor_var, "_sc"))),
  data = df_sf
)

# (B) Add topo variables
m_lm_topo <- lm(
  as.formula(paste(
    response_var, "~",
    paste0(predictor_var, "_sc"),
    "+ avg_lvt_sc + elvtn_r_sc + dmnnt_v"
  )),
  data = df_sf
)

# (C) Interaction model (tests whether topo modifies spectral signal)
m_inter <- lm(
  as.formula(paste(
    response_var, "~",
    paste0(predictor_var, "_sc"),
    "* dmnnt_v + avg_lvt_sc + elvtn_r_sc"
  )),
  data = df_sf
)

# (D) Nonlinear GAM (flexible relationship test)
m_gam <- gam(
  as.formula(paste(
    response_var, "~",
    paste0("s(", predictor_var, "_sc) +"),
    "s(avg_lvt_sc) + s(elvtn_r_sc) + dmnnt_v"
  )),
  data = df_sf,
  method = "REML"
)

# ---------------------------
# 5) MODEL COMPARISON
# ---------------------------

model_list <- list(
  lm_spectral = m_lm,
  lm_full     = m_lm_topo,
  interaction = m_inter,
  gam_model   = m_gam
)

aic_table <- AIC(m_lm, m_lm_topo, m_inter)

print(aic_table)

# ---------------------------
# 6) TEST: Does topo add explanatory power beyond spectral?
# ---------------------------

anova(m_lm, m_lm_topo, m_inter)

# partial F-tests
summary(m_lm_topo)

# ---------------------------
# 7) GAM diagnostics (nonlinearity test)
# ---------------------------

summary(m_gam)
plot(m_gam, pages = 1, shade = TRUE)

# ---------------------------
# 8) VISUAL: best model fit
# ---------------------------

df_sf$fitted_best <- fitted(m_inter)

ggplot(df_sf, aes_string(
  x = predictor_var,
  y = response_var
)) +
  geom_point(aes(color = dmnnt_v), alpha = 0.6) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(aes(y = fitted_best), color = "red", alpha = 0.4) +
  labs(
    title = "PD vs Spectral Heterogeneity",
    subtitle = "Black = linear fit | Red = topo-informed model",
    color = "Topo class"
  ) +
  theme_minimal()

# ---------------------------
# 9) DIAGNOSTIC: does topo explain residual variation?
# ---------------------------

df_sf$resid_lm <- resid(m_lm)

lm_topo_resid <- lm(resid_lm ~ avg_lvt_sc + elvtn_r_sc + dmnnt_v, data = df_sf)
summary(lm_topo_resid)

# ---------------------------
# 10) INTERPRETATION OUTPUT
# ---------------------------

cat("
============================================================
MODEL INTERPRETATION SUMMARY
============================================================

QUESTION:
Does topography explain additional variation in PD beyond spectral heterogeneity?

------------------------------------------------------------
1. BASE RELATIONSHIP
------------------------------------------------------------
Linear spectral model R²:
", summary(m_lm)$r.squared, "

------------------------------------------------------------
2. DOES TOPOGRAPHY ADD EXPLANATORY POWER?
------------------------------------------------------------
Full model R²:
", summary(m_lm_topo)$r.squared, "

ΔR²:
", summary(m_lm_topo)$r.squared - summary(m_lm)$r.squared, "

If ΔR² > 0.02 → topo meaningfully improves prediction
If ΔR² ~ 0     → spectral signal dominates

------------------------------------------------------------
3. NONLINEARITY CHECK (GAM)
------------------------------------------------------------
Look at:
- Effective degrees of freedom (edf)
- Significance of smooth terms

If edf >> 1 → strong nonlinear ecological structure
If edf ~ 1  → relationship is basically linear

------------------------------------------------------------
4. ECOLOGICAL INTERPRETATION
------------------------------------------------------------

- If lm_topo improves AIC:
  → Topography structures phylogenetic diversity beyond spectral heterogeneity

- If GAM is best:
  → PD responds nonlinearly to environmental gradients

- If only spectral is significant:
  → spectral heterogeneity likely proxies habitat complexity directly

============================================================
")
