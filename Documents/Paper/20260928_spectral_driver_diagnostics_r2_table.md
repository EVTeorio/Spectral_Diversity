# Table 7 Companion. R2 Values for Spectral Driver Diagnostics

R2 values for relationships between spectral heterogeneity metrics and illumination, brightness, elevation, and terrain drivers. This companion table mirrors Table 7 but reports explained variation rather than Pearson correlation coefficients.

## Table 7d. Regional Brightness R2 Values for Vector-Normalized PCA Metrics

| Spectral metric | Brightness region | 10 m R2 | 20 m R2 | 50 m R2 |
|---|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | Blue (450-495 nm) | 0.212*** | 0.060*** | 0.055* |
| Vector-normalized PCA alpha-hull area | Green (500-570 nm) | 0.233*** | 0.110*** | 0.222*** |
| Vector-normalized PCA alpha-hull area | Red (620-750 nm) | 0.016*** | 0.004 | 0.003 |
| Vector-normalized PCA alpha-hull area | Near infrared (750-998 nm) | 0.035*** | 0.134*** | 0.177*** |
| Vector-normalized PCA mean Euclidean distance | Blue (450-495 nm) | 0.081*** | 0.042*** | 0.019 |
| Vector-normalized PCA mean Euclidean distance | Green (500-570 nm) | 0.077*** | 0.089*** | 0.148*** |
| Vector-normalized PCA mean Euclidean distance | Red (620-750 nm) | 0.003* | 0.004 | 0.014 |
| Vector-normalized PCA mean Euclidean distance | Near infrared (750-998 nm) | 0.124*** | 0.127*** | 0.217*** |
| Vector-normalized PCA spectral Rao's Q | Blue (450-495 nm) | 0.077*** | 0.053*** | 0.007 |
| Vector-normalized PCA spectral Rao's Q | Green (500-570 nm) | 0.056*** | 0.074*** | 0.060* |
| Vector-normalized PCA spectral Rao's Q | Red (620-750 nm) | 0.001 | 0.002 | 0.013 |
| Vector-normalized PCA spectral Rao's Q | Near infrared (750-998 nm) | 0.043*** | 0.041*** | 0.151*** |

## Table 7e. Regional Brightness R2 Values for Raw PCA Metrics

| Spectral metric | Brightness region | 10 m R2 | 20 m R2 | 50 m R2 |
|---|---|---:|---:|---:|
| Raw PCA alpha-hull area | Blue (450-495 nm) | 0.454*** | 0.350*** | 0.275*** |
| Raw PCA alpha-hull area | Green (500-570 nm) | 0.575*** | 0.518*** | 0.557*** |
| Raw PCA alpha-hull area | Red (620-750 nm) | 0.328*** | 0.240*** | 0.121** |
| Raw PCA alpha-hull area | Near infrared (750-998 nm) | 0.085*** | 0.041*** | 0.001 |
| Raw PCA mean Euclidean distance | Blue (450-495 nm) | 0.293*** | 0.275*** | 0.211*** |
| Raw PCA mean Euclidean distance | Green (500-570 nm) | 0.443*** | 0.492*** | 0.533*** |
| Raw PCA mean Euclidean distance | Red (620-750 nm) | 0.248*** | 0.272*** | 0.157*** |
| Raw PCA mean Euclidean distance | Near infrared (750-998 nm) | 0.047*** | 0.065*** | 0.001 |
| Raw PCA spectral Rao's Q | Blue (450-495 nm) | 0.218*** | 0.213*** | 0.133** |
| Raw PCA spectral Rao's Q | Green (500-570 nm) | 0.284*** | 0.356*** | 0.368*** |
| Raw PCA spectral Rao's Q | Red (620-750 nm) | 0.135*** | 0.174*** | 0.081* |
| Raw PCA spectral Rao's Q | Near infrared (750-998 nm) | 0.009*** | 0.022** | 0.004 |

## Table 7f. Overall Brightness, Elevation, and Terrain R2 Values for Vector-Normalized PCA Metrics

| Spectral metric | Driver | 10 m R2 | 20 m R2 | 50 m R2 |
|---|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | 563 nm brightness | 0.189*** | 0.087*** | 0.186*** |
| Vector-normalized PCA alpha-hull area | Mean elevation | 0.032*** | 0.102*** | 0.236*** |
| Vector-normalized PCA alpha-hull area | TRI 5 x 5 | 0.008*** | 0.026*** | 0.058* |
| Vector-normalized PCA alpha-hull area | TRI 11 x 11 | 0.012*** | 0.033*** | 0.065* |
| Vector-normalized PCA alpha-hull area | TRI 21 x 21 | 0.018*** | 0.044*** | 0.079* |
| Vector-normalized PCA mean Euclidean distance | 563 nm brightness | 0.059*** | 0.075*** | 0.127** |
| Vector-normalized PCA mean Euclidean distance | Mean elevation | 0.041*** | 0.088*** | 0.232*** |
| Vector-normalized PCA mean Euclidean distance | TRI 5 x 5 | 0.005** | 0.016** | 0.016 |
| Vector-normalized PCA mean Euclidean distance | TRI 11 x 11 | 0.008*** | 0.020** | 0.022 |
| Vector-normalized PCA mean Euclidean distance | TRI 21 x 21 | 0.013*** | 0.028*** | 0.035 |
| Vector-normalized PCA spectral Rao's Q | 563 nm brightness | 0.046*** | 0.069*** | 0.054* |
| Vector-normalized PCA spectral Rao's Q | Mean elevation | 0.010*** | 0.026*** | 0.101** |
| Vector-normalized PCA spectral Rao's Q | TRI 5 x 5 | 0.001 | 0.003 | 0.001 |
| Vector-normalized PCA spectral Rao's Q | TRI 11 x 11 | 0.002 | 0.004 | 0.000 |
| Vector-normalized PCA spectral Rao's Q | TRI 21 x 21 | 0.003* | 0.006 | 0.001 |

Sources: `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv` and `reports/tables/final_research_direction/metric_driver_relationships.csv`.

Notes: Regional brightness values summarize retained-pixel brightness in blue (450-495 nm), green (500-570 nm), red (620-750 nm), and near-infrared (750-998 nm) regions. The 563 nm brightness variable is the mean retained-pixel reflectance at the shadow-mask wavelength. Complete-case sample sizes for vector-normalized PCA metrics were 1744 at 10 m, 436 at 20 m, and 74 at 50 m. Significance: *** p < .001; ** p < .01; * p < .05.
