
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(MASS)  # For loess and other non-linear models

# Define edge quads (to filter out later)
edge_quads <- c(
  "0", "1", "100", "2", "200", "3", "300", "4", "400", "5", "500", "6", "600",
  "7", "700", "8", "800", "9", "900", "10", "1000", "11", "1100", "12", "1200", "13",
  "14", "1300", "15", "1400", "16", "1500", "17", "1600", "18", "1700", "19", "1800", "20",
  "1900", "21", "22", "1901", "23", "24", "1903", "124", "1905", "224", "1906", "324", "1907",
  "424", "1908", "524", "1909", "624", "1910", "1911", "724", "1912", "824", "1913", "924", "1914",
  "1024", "1915", "1124", "1916", "1224", "1324", "1917", "1424", "1919", "1524", "1920", "1624", "1921",
  "1724", "1922", "1824", "1923", "1924", "1904", "1902", "1918"
)

# Load and clean data
setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")
df <- read.csv("Spectral_Diversity/Indices_SHPs/20m_spectral_sp.csv") %>%
  clean_names() %>%
  dplyr::filter(
    complete.cases(.),
    !name %in% edge_quads   # <-- comment this line to see in scatter plot.
  )

# Convert geometry to sf
df_sf <- st_as_sf(df, wkt = "geometry", crs = 4326)
pt <- st_point_on_surface(df_sf$geometry)
xy <- st_coordinates(pt)
df_sf <- df_sf %>%
  mutate(lon_mean = xy[, 1], lat_mean = xy[, 2])

#Options
spectral_vars <- c("global_pca_20m_masked_5nm",
                   "sa_entropy_20m_masked_5nm",
                   "sa_entropy_smooth",
                   "sa_entropy_smooth_masked",
                   "sa_entropy_smooth_masked_711",
                   "sum_band_entropy_20m_masked_5nm"
)

pd_vars <- c("afaith_pd", "faith_pd", "rao_pd", "rchnss_", "shnnn_d", "smpsn_d")

# Define predictor & response variables
response_var <- "faith_pd"
predictor_var <-  "sa_entropy_20m_masked_5nm"

# Exploratory plot
scatter_plot <- ggplot(df_sf, aes_string(x = predictor_var, y = response_var)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE) +  # Linear fit
  geom_smooth(method = "loess", color = "blue", se = FALSE) +    # LOESS (non-linear)
  labs(
    title = paste(response_var, "vs", predictor_var),
    subtitle = "Red = OLS regression | Blue = LOESS non-linear fit",
    x = predictor_var,
    y = response_var
  ) +
  theme_minimal()
print(scatter_plot)

# Step 1: Test Linearity with Residuals Plot
# Fit a linear model and get residuals
lm_model <- lm(rao_pd ~ sa_entropy_smooth_masked_711, data = df_sf)

# Residual plot for linearity check
ggplot(data = df_sf, aes(x = sa_entropy_smooth_masked_711, y = residuals(lm_model))) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Residuals Plot", x = predictor_var, y = "Residuals") +
  theme_minimal()

# Step 2: Test fit of different models
# Fit polynomial models of degree 2 and 3
poly2_model <- lm(rao_pd ~ poly(sa_entropy_smooth_masked_711, 2), data = df_sf)
poly3_model <- lm(rao_pd ~ poly(sa_entropy_smooth_masked_711, 3), data = df_sf)

# Compare models using AIC (Lower AIC is better)
AIC(lm_model, poly2_model, poly3_model)

# Step 3: Apply the best fitting model (Based on AIC or residuals)
best_model <- lm_model  # Default is linear model if no better fit
if (AIC(poly2_model) < AIC(best_model)) best_model <- poly2_model
if (AIC(poly3_model) < AIC(best_model)) best_model <- poly3_model

# Step 4: Visualize best fit model
ggplot(df_sf, aes_string(x = predictor_var, y = response_var)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "firebrick", se = FALSE) +  # Best fit (polynomial 2)
  labs(
    title = paste("Best Fit Model: ", response_var, "vs", predictor_var),
    subtitle = paste("Best Model: Polynomial Degree 2 (AIC =", round(AIC(poly2_model), 2), ")"),
    x = predictor_var,
    y = response_var
  ) +
  theme_minimal()

# Step 5: Check for statistical significance of the chosen model
summary(best_model)
