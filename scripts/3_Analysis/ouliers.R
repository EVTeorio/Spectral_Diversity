
library(dplyr)
library(tidyr)
library(vegan)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

edge_quads <- c(
  "0","1","100","2","200","3","300","4","400","5","500","6","600",
  "7","700","8","800","9","900","10","1000","11","1100","12","1200","13",
  "14","1300","15","1400","16","1500","17","1600","18","1700","19","1800","20",
  "1900","21","22","1901","23","24","1903","124","1905","224","1906","324","1907",
  "424","1908","524","1909","624","1910","1911","724","1912","824","1913","924","1914",
  "1024","1915","1124","1916","1224","1324","1917","1424","1919","1524","1920","1624","1921",
  "1724","1922","1824","1923","1924","1904","1902","1918"
)

df <- read.csv("Spectral_Diversity/Quad_Values/20m_spectral_sp.csv") %>%
  clean_names() %>%
  dplyr::filter(
    complete.cases(.),
    !name %in% edge_quads   # <-- comment this line to see in scatter plot.
  )

#-----------------------------
# 1. DEFINE VARIABLE GROUPS
#-----------------------------

spectral_vars <- c(
  "sa_entropy_smooth_masked_711",
  "sa_entropy_smooth_masked",
  "sa_entropy_smooth",
  "global_pca_20m_masked_5nm",
  "sa_entropy_20m_masked_5nm",
  "sum_band_entropy_20m_masked_5nm"
)

pd_vars <- c("afaith_pd", "faith_pd", "rao_pd")

biodiv_vars <- c(
  "afaith_pd", "faith_pd", "rao_pd", "rchnss",
  "shnnn_d", "smpsn_d", "evnnss"
)

all_vars <- c(spectral_vars, biodiv_vars)

#-----------------------------
# 2. OUTLIER FUNCTION
#-----------------------------

flag_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  
  return(x < lower | x > upper)
}

#-----------------------------
# SUMMARY TABLE WITH EDGE INFO
#-----------------------------
outlier_counts <- data.frame(
  variable = character(),
  n_outliers = numeric(),
  n_edge_outliers = numeric(),
  stringsAsFactors = FALSE
)

#-----------------------------
# LOOP
#-----------------------------
for (v in all_vars) {
  
  is_outlier <- flag_outliers(df[[v]])
  
  # Names of outliers
  outlier_names <- df$name[is_outlier]
  
  # Count total outliers
  n_out <- sum(is_outlier, na.rm = TRUE)
  
  # Count how many are in edge_quads
  n_edge <- sum(outlier_names %in% edge_quads)
  
  # Store
  outlier_counts <- rbind(
    outlier_counts,
    data.frame(
      variable = v,
      n_outliers = n_out,
      n_edge_outliers = n_edge
    )
  )
}

outlier_counts

# Topo Influence on spectral/ PD measures

results <- data.frame()

for (y in pd_vars) {
  for (x in spectral_vars) {
    
    formula <- as.formula(paste(y, "~", x, "* dmnnt_v"))
    
    model <- lm(formula, data = df)
    summ <- summary(model)
    
    # Extract interaction p-value(s)
    pvals <- coef(summ)[,4]
    interaction_terms <- grep(":", names(pvals), value = TRUE)
    
    if (length(interaction_terms) > 0) {
      p_int <- min(pvals[interaction_terms], na.rm = TRUE)
    } else {
      p_int <- NA
    }
    
    results <- rbind(results, data.frame(
      pd_variable = y,
      spectral_variable = x,
      interaction_p = p_int,
      r_squared = summ$r.squared
    ))
  }
}

results


library(ggplot2)

ggplot(df, aes(x = sa_entropy_20m_masked_5nm, y = afaith_pd, color = dmnnt_v)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

# Example for one PD variable
df$afaith_outlier <- flag_outliers(df$afaith_pd)

table(df$afaith_outlier, df$dmnnt_v)
fisher.test(table(df$afaith_outlier, df$dmnnt_v))

model <- lm(afaith_pd ~ sa_entropy_smooth * dmnnt_v, data = df)

# Extract residuals
df$resid <- residuals(model)

# Identify residual outliers
df$resid_outlier <- abs(scale(df$resid)) > 2

table(df$afaith_outlier, df$resid_outlier)
#####################################################

#######################################################

Y <- df %>%
  select(all_of(pd_vars)) %>%
  as.data.frame()

# Predictors
X <- df %>%
  select(
    sa_entropy_smooth_masked_711,
    sa_entropy_smooth_masked,
    sa_entropy_smooth,
    global_pca_20m_masked_5nm,
    sa_entropy_20m_masked_5nm,
    sum_band_entropy_20m_masked_5nm,
    avg_lvt,
    dmnnt_v
  ) %>%
  mutate(dmnnt_v = as.factor(dmnnt_v))

rda_model <- rda(Y ~ ., data = X)

summary(rda_model)
anova(rda_model)              # overall significance
anova(rda_model, by = "term") # variable importance

cca_model <- cca(Y ~ ., data = X)

summary(cca_model)
anova(cca_model)
anova(cca_model, by = "term")

plot(rda_model, scaling = 2)

library(ggplot2)

scores_sites <- scores(rda_model, display = "sites", scaling = 2)
scores_vars  <- scores(rda_model, display = "bp", scaling = 2)

sites_df <- as.data.frame(scores_sites)
sites_df$elevation <- df$avg_lvt
sites_df$topo <- df$dmnnt_v

vars_df <- as.data.frame(scores_vars)
vars_df$varnames <- rownames(vars_df)

ggplot() +
  geom_point(data = sites_df,
             aes(RDA1, RDA2, color = elevation),
             size = 3) +
  geom_segment(data = vars_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black") +
  geom_text(data = vars_df,
            aes(RDA1, RDA2, label = varnames),
            size = 3) +
  scale_color_viridis_c() +
  theme_minimal()
###############################
## Adding Outliers #########
###########################

plot(rda_model, scaling = 2)
points(rda_model, display = "sites", col = ifelse(df$afaith_outlier, "red", "black"))
df$rda1 <- scores(rda_model, display = "sites")[,1]

# Are outliers far from center?
t.test(df$rda1 ~ df$afaith_outlier)

site_scores <- scores(rda_model, display = "sites", scaling = 2)

df$rda1 <- site_scores[, 1]
flag_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  
  return(x < lower | x > upper)
}

df$rda1_outlier <- flag_outliers(df$rda1)

rda1_outlier_names <- df$name[df$rda1_outlier]
rda1_outlier_names
df$afaith_outlier

table(df$afaith_outlier, df$rda1_outlier)
pd_outliers <- df$name[df$afaith_outlier]

intersect(pd_outliers, rda1_outlier_names)
setdiff(pd_outliers, rda1_outlier_names)   # PD outliers NOT on RDA1 extremes
setdiff(rda1_outlier_names, pd_outliers)   # RDA1 extremes NOT PD outliers

###### Continuation #########################

#---------------------------------------------------
# 0. PREP
#---------------------------------------------------
df$dmnnt_v <- as.factor(df$dmnnt_v)
df$name <- as.character(df$name)

#---------------------------------------------------
# 1. EXTRACT RDA SCORES
#---------------------------------------------------
site_scores <- as.data.frame(scores(rda_model, display = "sites", scaling = 2))

df$rda1 <- site_scores[,1]
df$rda2 <- site_scores[,2]

#---------------------------------------------------
# 2. DEFINE OUTLIER FUNCTION (IQR RULE)
#---------------------------------------------------
flag_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  
  return(x < lower | x > upper)
}

# RDA1 outliers
df$rda1_outlier <- flag_outliers(df$rda1)

#---------------------------------------------------
# 3. ECOLOGICAL REGIMES (CLUSTERING)
#---------------------------------------------------
set.seed(123)
km <- kmeans(df[, c("rda1", "rda2")], centers = 4)
df$eco_cluster <- as.factor(km$cluster)

#---------------------------------------------------
# 4. VISUALIZATION: RDA SPACE
#---------------------------------------------------

# RDA colored by elevation
ggplot(df, aes(rda1, rda2)) +
  geom_point(aes(color = avg_lvt), size = 3) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "RDA space colored by elevation")

# RDA colored by topo
ggplot(df, aes(rda1, rda2)) +
  geom_point(aes(color = dmnnt_v), size = 3) +
  theme_minimal() +
  labs(title = "RDA space colored by topo type")

# RDA with outliers
df$outlier_group <- ifelse(df$afaith_outlier, "Outlier", "Normal")

ggplot(df, aes(rda1, rda2)) +
  geom_point(aes(color = outlier_group), size = 3) +
  theme_minimal() +
  labs(title = "PD outliers in RDA space")

#---------------------------------------------------
# 5. ECOLOGICAL REGIME VISUALIZATION
#---------------------------------------------------
ggplot(df, aes(rda1, rda2)) +
  geom_point(aes(color = eco_cluster), size = 3) +
  theme_minimal() +
  labs(title = "Ecological regimes (k-means on RDA space)")

#---------------------------------------------------
# 6. OUTLIERS vs REGIMES
#---------------------------------------------------
table(df$afaith_outlier, df$eco_cluster)
fisher.test(table(df$afaith_outlier, df$eco_cluster))

#---------------------------------------------------
# 7. ENVIRONMENT OF RDA1 EXTREMES
#---------------------------------------------------

# identify RDA1 extremes
df$rda1_outlier <- flag_outliers(df$rda1)

aggregate(avg_lvt ~ rda1_outlier, data = df, mean)
aggregate(avg_lvt ~ rda1_outlier, data = df, sd)

table(df$rda1_outlier, df$dmnnt_v)
prop.table(table(df$rda1_outlier, df$dmnnt_v), 1)

aggregate(sa_entropy_smooth_masked_711 ~ rda1_outlier, data = df, mean)
aggregate(sa_entropy_20m_masked_5nm ~ rda1_outlier, data = df, mean)

#---------------------------------------------------
# 8. CONNECT RDA1 EXTREMES TO PD OUTLIERS
#---------------------------------------------------
table(df$afaith_outlier, df$rda1_outlier)

pd_outliers <- df$name[df$afaith_outlier]
rda1_outliers <- df$name[df$rda1_outlier]

intersect(pd_outliers, rda1_outliers)

setdiff(pd_outliers, rda1_outliers)
setdiff(rda1_outliers, pd_outliers)

#---------------------------------------------------
# 9. OPTIONAL: BASIC SUMMARY STATS
#---------------------------------------------------
aggregate(rda1 ~ dmnnt_v, data = df, mean)
aggregate(rda1 ~ eco_cluster, data = df, mean)
aggregate(avg_lvt ~ eco_cluster, data = df, mean)

########################################################
## Plotting back to scatter plot 
######################################

#---------------------------------------------------
# USER SETTINGS
#---------------------------------------------------
response_var  <- "afaith_pd"
predictor_var <- "sa_entropy_20m_masked_5nm"

#---------------------------------------------------
# 1. BASE PREP (ensure factors/logicals exist)
#---------------------------------------------------
df$eco_cluster <- as.factor(df$eco_cluster)
df$dmnnt_v <- as.factor(df$dmnnt_v)
df$afaith_outlier <- as.logical(df$afaith_outlier)
df$rda1_outlier <- as.logical(df$rda1_outlier)

#---------------------------------------------------
# 2. BASE SCATTER FUNCTION (reusable)
#---------------------------------------------------
base_plot <- function(color_var, title_text) {
  ggplot(df, aes_string(x = predictor_var, y = response_var)) +
    geom_point(aes_string(color = color_var), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    theme_minimal() +
    labs(
      title = title_text,
      x = predictor_var,
      y = response_var
    )
}

#---------------------------------------------------
# 3. PLOT 1 — ECOLOGICAL CLUSTERS (REGIMES)
#---------------------------------------------------
p_cluster <- ggplot(df, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = eco_cluster), size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral Heterogeneity (colored by ecological regime)",
    color = "Eco cluster"
  )

print(p_cluster)

#---------------------------------------------------
# 4. PLOT 2 — ELEVATION GRADIENT
#---------------------------------------------------
p_elev <- ggplot(df, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = avg_lvt), size = 3, alpha = 0.85) +
  scale_color_viridis_c() +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral Heterogeneity (colored by elevation)",
    color = "Elevation"
  )

print(p_elev)

#---------------------------------------------------
# 5. PLOT 3 — TOPO TYPE STRUCTURE
#---------------------------------------------------
p_topo <- ggplot(df, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = dmnnt_v), size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral Heterogeneity (colored by topo type)",
    color = "Topo type"
  )

print(p_topo)

#---------------------------------------------------
# 6. PLOT 4 — OUTLIERS (FINAL STRUCTURE VIEW)
#---------------------------------------------------
p_outliers <- ggplot(df, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = afaith_outlier), size = 3, alpha = 0.85) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral Heterogeneity (highlighting outliers)",
    color = "PD outlier"
  )

print(p_outliers)

#---------------------------------------------------
# RDA1 OUTLIERS VIEW
#---------------------------------------------------
p_rda1_outliers <- ggplot(df, aes_string(x = predictor_var, y = response_var)) +
  geom_point(aes(color = rda1_outlier), size = 3, alpha = 0.85) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "blue")) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "PD vs Spectral Heterogeneity (RDA1 extremes)",
    color = "RDA1 outlier"
  )

print(p_rda1_outliers)

