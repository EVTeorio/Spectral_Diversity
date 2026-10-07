# Table 7. Spectral Driver Diagnostics for PCA Spectral Heterogeneity Metrics

Pearson correlations between spectral heterogeneity metrics and illumination, brightness, elevation, and terrain drivers. Values are shown for the primary vector-normalized PCA metrics and their raw PCA counterparts where useful for interpreting brightness dominance.

## Table 7a. Regional Brightness Correlations for Vector-Normalized PCA Metrics

| Spectral metric | Brightness region | 10 m r | 20 m r | 50 m r |
|---|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | Blue (450-495 nm) | 0.461*** | 0.245*** | 0.235* |
| Vector-normalized PCA alpha-hull area | Green (500-570 nm) | 0.483*** | 0.332*** | 0.471*** |
| Vector-normalized PCA alpha-hull area | Red (620-750 nm) | 0.126*** | -0.066 | -0.058 |
| Vector-normalized PCA alpha-hull area | Near infrared (750-998 nm) | -0.187*** | -0.366*** | -0.421*** |
| Vector-normalized PCA mean Euclidean distance | Blue (450-495 nm) | 0.285*** | 0.205*** | 0.137 |
| Vector-normalized PCA mean Euclidean distance | Green (500-570 nm) | 0.278*** | 0.298*** | 0.385*** |
| Vector-normalized PCA mean Euclidean distance | Red (620-750 nm) | -0.054* | -0.060 | -0.120 |
| Vector-normalized PCA mean Euclidean distance | Near infrared (750-998 nm) | -0.352*** | -0.357*** | -0.466*** |
| Vector-normalized PCA spectral Rao's Q | Blue (450-495 nm) | 0.278*** | 0.230*** | 0.081 |
| Vector-normalized PCA spectral Rao's Q | Green (500-570 nm) | 0.238*** | 0.272*** | 0.245* |
| Vector-normalized PCA spectral Rao's Q | Red (620-750 nm) | 0.032 | 0.045 | -0.114 |
| Vector-normalized PCA spectral Rao's Q | Near infrared (750-998 nm) | -0.208*** | -0.203*** | -0.389*** |

## Table 7b. Regional Brightness Correlations for Raw PCA Metrics

| Spectral metric | Brightness region | 10 m r | 20 m r | 50 m r |
|---|---|---:|---:|---:|
| Raw PCA alpha-hull area | Blue (450-495 nm) | 0.674*** | 0.592*** | 0.524*** |
| Raw PCA alpha-hull area | Green (500-570 nm) | 0.758*** | 0.720*** | 0.746*** |
| Raw PCA alpha-hull area | Red (620-750 nm) | 0.573*** | 0.489*** | 0.348** |
| Raw PCA alpha-hull area | Near infrared (750-998 nm) | 0.292*** | 0.202*** | -0.022 |
| Raw PCA mean Euclidean distance | Blue (450-495 nm) | 0.542*** | 0.524*** | 0.459*** |
| Raw PCA mean Euclidean distance | Green (500-570 nm) | 0.666*** | 0.702*** | 0.730*** |
| Raw PCA mean Euclidean distance | Red (620-750 nm) | 0.498*** | 0.522*** | 0.396*** |
| Raw PCA mean Euclidean distance | Near infrared (750-998 nm) | 0.217*** | 0.255*** | 0.036 |
| Raw PCA spectral Rao's Q | Blue (450-495 nm) | 0.467*** | 0.462*** | 0.365** |
| Raw PCA spectral Rao's Q | Green (500-570 nm) | 0.533*** | 0.597*** | 0.606*** |
| Raw PCA spectral Rao's Q | Red (620-750 nm) | 0.367*** | 0.418*** | 0.284* |
| Raw PCA spectral Rao's Q | Near infrared (750-998 nm) | 0.095*** | 0.148** | -0.063 |

## Table 7c. Overall Brightness, Elevation, and Terrain Correlations for Vector-Normalized PCA Metrics

| Spectral metric | Driver | 10 m r | 20 m r | 50 m r |
|---|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | 563 nm brightness | 0.435*** | 0.296*** | 0.432*** |
| Vector-normalized PCA alpha-hull area | Mean elevation | -0.178*** | -0.319*** | -0.485*** |
| Vector-normalized PCA alpha-hull area | TRI 5 x 5 | -0.091*** | -0.161*** | -0.242* |
| Vector-normalized PCA alpha-hull area | TRI 11 x 11 | -0.111*** | -0.182*** | -0.255* |
| Vector-normalized PCA alpha-hull area | TRI 21 x 21 | -0.135*** | -0.209*** | -0.281* |
| Vector-normalized PCA mean Euclidean distance | 563 nm brightness | 0.242*** | 0.274*** | 0.357** |
| Vector-normalized PCA mean Euclidean distance | Mean elevation | -0.202*** | -0.296*** | -0.482*** |
| Vector-normalized PCA mean Euclidean distance | TRI 5 x 5 | -0.071** | -0.125** | -0.127 |
| Vector-normalized PCA mean Euclidean distance | TRI 11 x 11 | -0.089*** | -0.142** | -0.148 |
| Vector-normalized PCA mean Euclidean distance | TRI 21 x 21 | -0.112*** | -0.166*** | -0.187 |
| Vector-normalized PCA spectral Rao's Q | 563 nm brightness | 0.213*** | 0.262*** | 0.233* |
| Vector-normalized PCA spectral Rao's Q | Mean elevation | -0.098*** | -0.160*** | -0.318** |
| Vector-normalized PCA spectral Rao's Q | TRI 5 x 5 | -0.027 | -0.053 | 0.035 |
| Vector-normalized PCA spectral Rao's Q | TRI 11 x 11 | -0.039 | -0.065 | 0.009 |
| Vector-normalized PCA spectral Rao's Q | TRI 21 x 21 | -0.052* | -0.080 | -0.036 |

Sources: `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv` and `reports/tables/final_research_direction/metric_driver_relationships.csv`.

Notes: Regional brightness values summarize retained-pixel brightness in blue (450-495 nm), green (500-570 nm), red (620-750 nm), and near-infrared (750-998 nm) regions. The 563 nm brightness variable is the mean retained-pixel reflectance at the shadow-mask wavelength. Complete-case sample sizes for vector-normalized PCA metrics were 1744 at 10 m, 436 at 20 m, and 74 at 50 m. Significance: *** p < .001; ** p < .01; * p < .05.
