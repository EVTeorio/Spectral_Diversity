# Edge Quadrat And Bootstrap Sensitivity For SV-Diversity Correlations

Date: 2026-07-25

## Purpose

This analysis compares how edge quadrats affect correlations between diversity measures and spectral variation at the 10 m and 20 m scales. It also checks whether bootstrap variation in spectral angle entropy changes the apparent strength of biodiversity-SV correlations.

## Inputs

- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv`

## Measures

- Diversity measures: `phy_faith`, `phy_rao`, `phy_afaith`, and `sp_shannon`.
- Spectral variation measures: `spec_sa` and all standardized PCA spectral variation fields in the current analysis dataset.

## Edge Quadrat Comparison

Each pairwise correlation was run three ways within 10 m and 20 m quadrats: all quadrats, non-edge quadrats only, and edge quadrats only. Scatter plots show all quadrats, highlight edge quadrats, and include all-quadrat, non-edge, and edge-only regression lines.

### Primary Measure Summary

| Scale | SV measure | Diversity measure | r all | R2 all | r non-edge | R2 non-edge | r edge-only | R2 edge-only | delta r | edge-only minus non-edge r |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | SA entropy | Faith's PD | 0.098 | 0.010 | 0.088 | 0.008 | 0.136 | 0.018 | -0.011 | 0.048 |
| 10m | SA entropy | Phylogenetic Rao's Q | 0.115 | 0.013 | 0.106 | 0.011 | 0.148 | 0.022 | -0.009 | 0.042 |
| 10m | SA entropy | Abundance-weighted Faith's PD | 0.106 | 0.011 | 0.102 | 0.010 | 0.111 | 0.012 | -0.004 | 0.009 |
| 10m | SA entropy | Shannon diversity | -0.027 | 0.001 | -0.051 | 0.003 | 0.116 | 0.014 | -0.023 | 0.167 |
| 10m | Std PCA mean distance | Faith's PD | 0.119 | 0.014 | 0.098 | 0.010 | 0.212 | 0.045 | -0.021 | 0.114 |
| 10m | Std PCA mean distance | Phylogenetic Rao's Q | 0.140 | 0.020 | 0.121 | 0.015 | 0.232 | 0.054 | -0.019 | 0.111 |
| 10m | Std PCA mean distance | Abundance-weighted Faith's PD | 0.122 | 0.015 | 0.106 | 0.011 | 0.215 | 0.046 | -0.016 | 0.109 |
| 10m | Std PCA mean distance | Shannon diversity | -0.032 | 0.001 | -0.049 | 0.002 | 0.068 | 0.005 | -0.017 | 0.118 |
| 20m | SA entropy | Faith's PD | 0.170 | 0.029 | 0.140 | 0.020 | 0.329 | 0.108 | -0.029 | 0.188 |
| 20m | SA entropy | Phylogenetic Rao's Q | 0.203 | 0.041 | 0.180 | 0.033 | 0.334 | 0.112 | -0.022 | 0.154 |
| 20m | SA entropy | Abundance-weighted Faith's PD | 0.181 | 0.033 | 0.166 | 0.027 | 0.250 | 0.063 | -0.016 | 0.085 |
| 20m | SA entropy | Shannon diversity | -0.021 | 0.000 | -0.065 | 0.004 | 0.326 | 0.106 | -0.044 | 0.391 |
| 20m | Std PCA mean distance | Faith's PD | 0.240 | 0.058 | 0.222 | 0.049 | 0.265 | 0.070 | -0.018 | 0.043 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q | 0.274 | 0.075 | 0.249 | 0.062 | 0.382 | 0.146 | -0.025 | 0.133 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD | 0.301 | 0.091 | 0.291 | 0.085 | 0.281 | 0.079 | -0.010 | -0.010 |
| 20m | Std PCA mean distance | Shannon diversity | -0.036 | 0.001 | -0.075 | 0.006 | 0.200 | 0.040 | -0.039 | 0.275 |

### Largest Edge-Removal Changes

| Scale | SV measure | Diversity measure | n all | edge n all | r all | r non-edge | delta r |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20m | Std PCA median distance | Shannon diversity |  436 |  76 | -0.011 | -0.056 | -0.045 |
| 20m | SA entropy | Shannon diversity |  485 |  80 | -0.021 | -0.065 | -0.044 |
| 20m | Std PCA mean distance | Shannon diversity |  436 |  76 | -0.036 | -0.075 | -0.039 |
| 20m | Std PCA alpha-hull area | Shannon diversity |  436 |  76 | -0.024 | -0.061 | -0.036 |
| 20m | Std PCA Rao's Q | Phylogenetic Rao's Q |  436 |  76 | 0.132 | 0.102 | -0.030 |
| 20m | Std PCA Rao's Q | Shannon diversity |  436 |  76 | -0.061 | -0.091 | -0.030 |
| 20m | SA entropy | Faith's PD |  485 |  80 | 0.170 | 0.140 | -0.029 |
| 10m | Std PCA convex-hull area | Abundance-weighted Faith's PD | 1744 | 304 | -0.004 | -0.030 | -0.026 |
| 10m | Std PCA 3D hull volume | Abundance-weighted Faith's PD | 1744 | 304 | -0.003 | -0.029 | -0.026 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q |  436 |  76 | 0.274 | 0.249 | -0.025 |

### Strongest Edge-Only Correlations

| Scale | SV measure | Diversity measure | n edge-only | r edge-only | R2 edge-only | F p-value edge-only | r non-edge |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20m | Std PCA median distance | Phylogenetic Rao's Q | 76 | 0.407 | 0.165 | 2.66e-04 | 0.327 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q | 76 | 0.382 | 0.146 | 6.66e-04 | 0.249 |
| 20m | SA entropy | Phylogenetic Rao's Q | 80 | 0.334 | 0.112 | 0.002 | 0.180 |
| 20m | SA entropy | Faith's PD | 80 | 0.329 | 0.108 | 0.003 | 0.140 |
| 20m | SA entropy | Shannon diversity | 80 | 0.326 | 0.106 | 0.003 | -0.065 |
| 20m | Std PCA alpha-hull area | Phylogenetic Rao's Q | 76 | 0.307 | 0.094 | 0.007 | 0.325 |
| 20m | Std PCA Rao's Q | Phylogenetic Rao's Q | 76 | 0.302 | 0.091 | 0.008 | 0.102 |
| 20m | Std PCA median distance | Faith's PD | 76 | 0.302 | 0.091 | 0.008 | 0.277 |
| 20m | Std PCA median distance | Abundance-weighted Faith's PD | 76 | 0.301 | 0.091 | 0.008 | 0.379 |
| 20m | Std PCA alpha-hull area | Abundance-weighted Faith's PD | 76 | 0.292 | 0.085 | 0.010 | 0.381 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD | 76 | 0.281 | 0.079 | 0.014 | 0.291 |
| 20m | Std PCA median distance | Shannon diversity | 76 | 0.271 | 0.073 | 0.018 | -0.056 |

## Edge Versus Non-Edge Range Check

This check compares the spread of diversity and spectral variables in edge and non-edge quadrats. Smaller edge-only ranges or standard deviations can make the scatterplots look less diffuse even when the underlying relationship is not fundamentally stronger.

| Scale | Variable | n edge | edge range | non-edge range | edge/non-edge range | edge SD | non-edge SD | edge/non-edge SD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | Faith's PD | 344 | 1300.666 | 1282.603 | 1.014 | 205.204 | 228.109 | 0.900 |
| 10m | Phylogenetic Rao's Q | 344 | 338.393 | 373.607 | 0.906 | 60.117 | 71.676 | 0.839 |
| 10m | Abundance-weighted Faith's PD | 344 | 819.200 | 1312.099 | 0.624 | 140.370 | 186.418 | 0.753 |
| 10m | Shannon diversity | 344 | 1.873 | 2.031 | 0.922 | 0.396 | 0.410 | 0.964 |
| 10m | SA entropy | 316 | 0.908 | 1.273 | 0.714 | 0.095 | 0.111 | 0.854 |
| 10m | Std PCA mean distance | 304 | 18.663 | 27.167 | 0.687 | 1.480 | 1.876 | 0.789 |
| 20m | Faith's PD |  86 | 1251.011 | 1684.947 | 0.742 | 239.832 | 286.908 | 0.836 |
| 20m | Phylogenetic Rao's Q |  86 | 217.732 | 379.908 | 0.573 | 43.117 | 58.391 | 0.738 |
| 20m | Abundance-weighted Faith's PD |  86 | 1624.601 | 3451.101 | 0.471 | 292.539 | 436.330 | 0.670 |
| 20m | Shannon diversity |  86 | 1.822 | 2.401 | 0.759 | 0.352 | 0.388 | 0.908 |
| 20m | SA entropy |  80 | 0.397 | 0.941 | 0.422 | 0.072 | 0.102 | 0.709 |
| 20m | Std PCA mean distance |  76 | 6.649 | 11.313 | 0.588 | 1.013 | 1.439 | 0.704 |

## Equal-Sample-Size Non-Edge Resampling

For each SV-diversity pair, non-edge quadrats were randomly downsampled to the matching edge-only sample size for 2000 iterations. The edge-only correlation was then compared with the distribution of same-sized non-edge correlations.

| Scale | SV measure | Diversity measure | n edge | edge r | non-edge equal-n mean r | non-edge equal-n 2.5% r | non-edge equal-n 97.5% r | abs-tail prob | same-direction tail prob |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | SA entropy | Faith's PD | 316 | 0.136 | 0.089 | -0.005 | 0.186 | 0.166 | 0.166 |
| 10m | SA entropy | Phylogenetic Rao's Q | 316 | 0.148 | 0.106 | 0.013 | 0.204 | 0.203 | 0.203 |
| 10m | SA entropy | Abundance-weighted Faith's PD | 316 | 0.111 | 0.102 | -0.007 | 0.212 | 0.419 | 0.419 |
| 10m | SA entropy | Shannon diversity | 316 | 0.116 | -0.051 | -0.144 | 0.043 | 0.081 | 0.000 |
| 10m | Std PCA mean distance | Faith's PD | 304 | 0.212 | 0.106 | -0.009 | 0.236 | 0.048 | 0.048 |
| 10m | Std PCA mean distance | Phylogenetic Rao's Q | 304 | 0.232 | 0.127 | 0.015 | 0.253 | 0.046 | 0.046 |
| 10m | Std PCA mean distance | Abundance-weighted Faith's PD | 304 | 0.215 | 0.113 | -0.005 | 0.247 | 0.064 | 0.064 |
| 10m | Std PCA mean distance | Shannon diversity | 304 | 0.068 | -0.048 | -0.147 | 0.048 | 0.353 | 0.012 |
| 20m | SA entropy | Faith's PD |  80 | 0.329 | 0.147 | -0.047 | 0.340 | 0.033 | 0.033 |
| 20m | SA entropy | Phylogenetic Rao's Q |  80 | 0.334 | 0.182 | -0.018 | 0.382 | 0.073 | 0.073 |
| 20m | SA entropy | Abundance-weighted Faith's PD |  80 | 0.250 | 0.166 | -0.043 | 0.357 | 0.220 | 0.220 |
| 20m | SA entropy | Shannon diversity |  80 | 0.326 | -0.064 | -0.292 | 0.164 | 0.011 | 0.000 |
| 20m | Std PCA mean distance | Faith's PD |  76 | 0.265 | 0.229 | 0.007 | 0.446 | 0.389 | 0.389 |
| 20m | Std PCA mean distance | Phylogenetic Rao's Q |  76 | 0.382 | 0.257 | 0.046 | 0.469 | 0.128 | 0.128 |
| 20m | Std PCA mean distance | Abundance-weighted Faith's PD |  76 | 0.281 | 0.299 | 0.085 | 0.513 | 0.574 | 0.574 |
| 20m | Std PCA mean distance | Shannon diversity |  76 | 0.200 | -0.074 | -0.273 | 0.123 | 0.115 | 0.003 |

## Environmental Context

This section checks whether elevation or topographic roughness differs between edge and non-edge quadrats, correlates with spectral variation, or reduces the diversity-SV relationship when included in the same model. Environmental adjusted models use `env_elev`, `env_tri5`, `env_tri11`, and `env_tri21` together.

### Edge Versus Non-Edge Environmental Distributions

| Scale | Variable | edge mean | non-edge mean | edge minus non-edge | edge/non-edge SD |
| --- | --- | --- | --- | --- | --- |
| 10m | Elevation | 230.859 | 223.700 | 7.160 | 1.485 |
| 10m | TRI 5x5 | 2.021 | 1.985 | 0.036 | 0.900 |
| 10m | TRI 11x11 | 9.416 | 9.136 | 0.281 | 0.905 |
| 10m | TRI 21x21 | 33.229 | 31.757 | 1.472 | 0.942 |
| 20m | Elevation | 230.859 | 223.700 | 7.160 | 1.493 |
| 20m | TRI 5x5 | 2.021 | 1.985 | 0.036 | 0.954 |
| 20m | TRI 11x11 | 9.416 | 9.136 | 0.281 | 0.953 |
| 20m | TRI 21x21 | 33.229 | 31.757 | 1.472 | 0.976 |

### Strongest Environment-SV Correlations

| Scale | Group | Environment | Variable | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20m | non_edge | Elevation | Std PCA mean distance |  360 | -0.306 | 0.094 | 3.07e-09 |
| 20m | all_quads | Elevation | Std PCA mean distance |  436 | -0.296 | 0.088 | 2.93e-10 |
| 20m | non_edge | Elevation | SA entropy |  405 | -0.294 | 0.087 | 1.53e-09 |
| 20m | all_quads | Elevation | SA entropy |  485 | -0.288 | 0.083 | 9.90e-11 |
| 20m | edge_only | Elevation | SA entropy |   80 | -0.281 | 0.079 | 0.012 |
| 20m | edge_only | Elevation | Std PCA mean distance |   76 | -0.258 | 0.067 | 0.024 |
| 10m | non_edge | Elevation | Std PCA mean distance | 1440 | -0.209 | 0.044 | 1.06e-15 |
| 10m | edge_only | Elevation | SA entropy |  316 | -0.204 | 0.042 | 2.59e-04 |
| 10m | all_quads | Elevation | Std PCA mean distance | 1744 | -0.202 | 0.041 | 1.82e-17 |
| 10m | all_quads | Elevation | SA entropy | 1903 | -0.183 | 0.034 | 8.03e-16 |
| 10m | non_edge | Elevation | SA entropy | 1587 | -0.177 | 0.031 | 1.13e-12 |
| 10m | edge_only | Elevation | Std PCA mean distance |  304 | -0.172 | 0.029 | 0.003 |
| 20m | all_quads | TRI 21x21 | Std PCA mean distance |  436 | -0.166 | 0.028 | 4.89e-04 |
| 20m | non_edge | TRI 21x21 | Std PCA mean distance |  360 | -0.165 | 0.027 | 0.002 |
| 20m | edge_only | TRI 21x21 | Std PCA mean distance |   76 | -0.159 | 0.025 | 0.170 |
| 20m | edge_only | TRI 21x21 | SA entropy |   80 | -0.159 | 0.025 | 0.159 |

### Diversity Models Before And After Environmental Adjustment

| Scale | Group | SV measure | Diversity measure | n | simple R2 | env-only R2 | adjusted R2 | diversity incremental R2 | diversity p adjusted |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | non_edge | SA entropy | Faith's PD | 1587 | 0.008 | 0.032 | 0.037 | 0.005 | 0.004 |
| 10m | non_edge | SA entropy | Phylogenetic Rao's Q | 1587 | 0.011 | 0.032 | 0.041 | 0.009 | 1.38e-04 |
| 10m | non_edge | SA entropy | Abundance-weighted Faith's PD | 1587 | 0.010 | 0.032 | 0.040 | 0.008 | 2.71e-04 |
| 10m | non_edge | SA entropy | Shannon diversity | 1587 | 0.003 | 0.032 | 0.034 | 0.001 | 0.156 |
| 10m | non_edge | Std PCA mean distance | Faith's PD | 1440 | 0.010 | 0.056 | 0.061 | 0.005 | 0.007 |
| 10m | non_edge | Std PCA mean distance | Phylogenetic Rao's Q | 1440 | 0.015 | 0.056 | 0.065 | 0.009 | 2.47e-04 |
| 10m | non_edge | Std PCA mean distance | Abundance-weighted Faith's PD | 1440 | 0.011 | 0.056 | 0.062 | 0.006 | 0.003 |
| 10m | non_edge | Std PCA mean distance | Shannon diversity | 1440 | 0.002 | 0.056 | 0.057 | 0.001 | 0.292 |
| 10m | edge_only | SA entropy | Faith's PD |  316 | 0.018 | 0.045 | 0.054 | 0.009 | 0.095 |
| 10m | edge_only | SA entropy | Phylogenetic Rao's Q |  316 | 0.022 | 0.045 | 0.059 | 0.014 | 0.035 |
| 10m | edge_only | SA entropy | Abundance-weighted Faith's PD |  316 | 0.012 | 0.045 | 0.054 | 0.009 | 0.095 |
| 10m | edge_only | SA entropy | Shannon diversity |  316 | 0.014 | 0.045 | 0.055 | 0.010 | 0.078 |
| 10m | edge_only | Std PCA mean distance | Faith's PD |  304 | 0.045 | 0.040 | 0.080 | 0.040 | 3.56e-04 |
| 10m | edge_only | Std PCA mean distance | Phylogenetic Rao's Q |  304 | 0.054 | 0.040 | 0.089 | 0.049 | 7.90e-05 |
| 10m | edge_only | Std PCA mean distance | Abundance-weighted Faith's PD |  304 | 0.046 | 0.040 | 0.088 | 0.048 | 8.87e-05 |
| 10m | edge_only | Std PCA mean distance | Shannon diversity |  304 | 0.005 | 0.040 | 0.045 | 0.005 | 0.221 |
| 20m | non_edge | SA entropy | Faith's PD |  405 | 0.020 | 0.090 | 0.094 | 0.004 | 0.176 |
| 20m | non_edge | SA entropy | Phylogenetic Rao's Q |  405 | 0.033 | 0.090 | 0.105 | 0.016 | 0.008 |
| 20m | non_edge | SA entropy | Abundance-weighted Faith's PD |  405 | 0.027 | 0.090 | 0.102 | 0.012 | 0.022 |
| 20m | non_edge | SA entropy | Shannon diversity |  405 | 0.004 | 0.090 | 0.092 | 0.002 | 0.302 |
| 20m | non_edge | Std PCA mean distance | Faith's PD |  360 | 0.049 | 0.117 | 0.132 | 0.015 | 0.015 |
| 20m | non_edge | Std PCA mean distance | Phylogenetic Rao's Q |  360 | 0.062 | 0.117 | 0.142 | 0.025 | 0.001 |
| 20m | non_edge | Std PCA mean distance | Abundance-weighted Faith's PD |  360 | 0.085 | 0.117 | 0.161 | 0.044 | 2.23e-05 |
| 20m | non_edge | Std PCA mean distance | Shannon diversity |  360 | 0.006 | 0.117 | 0.121 | 0.004 | 0.229 |
| 20m | edge_only | SA entropy | Faith's PD |   80 | 0.108 | 0.119 | 0.160 | 0.041 | 0.060 |
| 20m | edge_only | SA entropy | Phylogenetic Rao's Q |   80 | 0.112 | 0.119 | 0.172 | 0.052 | 0.034 |
| 20m | edge_only | SA entropy | Abundance-weighted Faith's PD |   80 | 0.063 | 0.119 | 0.145 | 0.026 | 0.139 |
| 20m | edge_only | SA entropy | Shannon diversity |   80 | 0.106 | 0.119 | 0.164 | 0.044 | 0.051 |
| 20m | edge_only | Std PCA mean distance | Faith's PD |   76 | 0.070 | 0.068 | 0.119 | 0.051 | 0.048 |
| 20m | edge_only | Std PCA mean distance | Phylogenetic Rao's Q |   76 | 0.146 | 0.068 | 0.182 | 0.114 | 0.003 |
| 20m | edge_only | Std PCA mean distance | Abundance-weighted Faith's PD |   76 | 0.079 | 0.068 | 0.136 | 0.067 | 0.022 |
| 20m | edge_only | Std PCA mean distance | Shannon diversity |   76 | 0.040 | 0.068 | 0.089 | 0.021 | 0.207 |

## Bootstrap Variation Sensitivity

Bootstrap sensitivity was evaluated for `spec_sa`, because the bootstrap diagnostics table describes uncertainty in spectral angle entropy. The analysis compares non-edge correlations across all bootstrap-valid quadrats, low-CV subsets, high-CV subsets, and scale-specific bootstrap-SD tertiles.

| Scale | Subset | Diversity measure | n | mean boot SD | mean boot CV | r | R2 | F p-value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | non_edge_all_bootstrap | Faith's PD | 1554 | 0.0668 | 0.0236 | 0.099 | 0.010 | 9.55e-05 |
| 10m | boot_cv_le_5pct | Faith's PD | 1254 | 0.0241 | 0.0086 | 0.120 | 0.014 | 2.20e-05 |
| 10m | boot_cv_gt_5pct | Faith's PD |  300 | 0.2453 | 0.0863 | 0.016 | 0.000 | 0.785 |
| 10m | low_boot_sd | Faith's PD |  528 | 0.0082 | 0.0029 | 0.174 | 0.030 | 5.71e-05 |
| 10m | high_boot_sd | Faith's PD |  512 | 0.1839 | 0.0649 | 0.012 | 0.000 | 0.789 |
| 10m | non_edge_all_bootstrap | Phylogenetic Rao's Q | 1554 | 0.0668 | 0.0236 | 0.118 | 0.014 | 3.02e-06 |
| 10m | boot_cv_le_5pct | Phylogenetic Rao's Q | 1254 | 0.0241 | 0.0086 | 0.147 | 0.022 | 1.71e-07 |
| 10m | boot_cv_gt_5pct | Phylogenetic Rao's Q |  300 | 0.2453 | 0.0863 | 0.009 | 0.000 | 0.874 |
| 10m | low_boot_sd | Phylogenetic Rao's Q |  528 | 0.0082 | 0.0029 | 0.207 | 0.043 | 1.67e-06 |
| 10m | high_boot_sd | Phylogenetic Rao's Q |  512 | 0.1839 | 0.0649 | 0.001 | 0.000 | 0.977 |
| 10m | non_edge_all_bootstrap | Abundance-weighted Faith's PD | 1554 | 0.0668 | 0.0236 | 0.116 | 0.013 | 4.86e-06 |
| 10m | boot_cv_le_5pct | Abundance-weighted Faith's PD | 1254 | 0.0241 | 0.0086 | 0.137 | 0.019 | 1.03e-06 |
| 10m | boot_cv_gt_5pct | Abundance-weighted Faith's PD |  300 | 0.2453 | 0.0863 | 0.035 | 0.001 | 0.541 |
| 10m | low_boot_sd | Abundance-weighted Faith's PD |  528 | 0.0082 | 0.0029 | 0.142 | 0.020 | 0.001 |
| 10m | high_boot_sd | Abundance-weighted Faith's PD |  512 | 0.1839 | 0.0649 | 0.009 | 0.000 | 0.834 |
| 10m | non_edge_all_bootstrap | Shannon diversity | 1554 | 0.0668 | 0.0236 | -0.041 | 0.002 | 0.105 |
| 10m | boot_cv_le_5pct | Shannon diversity | 1254 | 0.0241 | 0.0086 | -0.041 | 0.002 | 0.144 |
| 10m | boot_cv_gt_5pct | Shannon diversity |  300 | 0.2453 | 0.0863 | -0.042 | 0.002 | 0.471 |
| 10m | low_boot_sd | Shannon diversity |  528 | 0.0082 | 0.0029 | -0.008 | 0.000 | 0.854 |
| 10m | high_boot_sd | Shannon diversity |  512 | 0.1839 | 0.0649 | -0.069 | 0.005 | 0.119 |
| 20m | non_edge_all_bootstrap | Faith's PD |  398 | 0.0521 | 0.0185 | 0.146 | 0.021 | 0.003 |
| 20m | boot_cv_le_5pct | Faith's PD |  349 | 0.0265 | 0.0093 | 0.149 | 0.022 | 0.005 |
| 20m | boot_cv_gt_5pct | Faith's PD |   49 | 0.2347 | 0.0837 | 0.037 | 0.001 | 0.803 |
| 20m | low_boot_sd | Faith's PD |  134 | 0.0107 | 0.0037 | 0.040 | 0.002 | 0.643 |
| 20m | high_boot_sd | Faith's PD |  129 | 0.1363 | 0.0484 | 0.057 | 0.003 | 0.518 |
| 20m | non_edge_all_bootstrap | Phylogenetic Rao's Q |  398 | 0.0521 | 0.0185 | 0.194 | 0.038 | 9.70e-05 |
| 20m | boot_cv_le_5pct | Phylogenetic Rao's Q |  349 | 0.0265 | 0.0093 | 0.219 | 0.048 | 3.65e-05 |
| 20m | boot_cv_gt_5pct | Phylogenetic Rao's Q |   49 | 0.2347 | 0.0837 | -0.058 | 0.003 | 0.690 |
| 20m | low_boot_sd | Phylogenetic Rao's Q |  134 | 0.0107 | 0.0037 | 0.154 | 0.024 | 0.077 |
| 20m | high_boot_sd | Phylogenetic Rao's Q |  129 | 0.1363 | 0.0484 | 0.068 | 0.005 | 0.443 |
| 20m | non_edge_all_bootstrap | Abundance-weighted Faith's PD |  398 | 0.0521 | 0.0185 | 0.179 | 0.032 | 3.37e-04 |
| 20m | boot_cv_le_5pct | Abundance-weighted Faith's PD |  349 | 0.0265 | 0.0093 | 0.202 | 0.041 | 1.40e-04 |
| 20m | boot_cv_gt_5pct | Abundance-weighted Faith's PD |   49 | 0.2347 | 0.0837 | -0.100 | 0.010 | 0.495 |
| 20m | low_boot_sd | Abundance-weighted Faith's PD |  134 | 0.0107 | 0.0037 | 0.168 | 0.028 | 0.052 |
| 20m | high_boot_sd | Abundance-weighted Faith's PD |  129 | 0.1363 | 0.0484 | 0.007 | 0.000 | 0.933 |
| 20m | non_edge_all_bootstrap | Shannon diversity |  398 | 0.0521 | 0.0185 | -0.071 | 0.005 | 0.157 |
| 20m | boot_cv_le_5pct | Shannon diversity |  349 | 0.0265 | 0.0093 | -0.047 | 0.002 | 0.376 |
| 20m | boot_cv_gt_5pct | Shannon diversity |   49 | 0.2347 | 0.0837 | -0.225 | 0.050 | 0.121 |
| 20m | low_boot_sd | Shannon diversity |  134 | 0.0107 | 0.0037 | -0.069 | 0.005 | 0.429 |
| 20m | high_boot_sd | Shannon diversity |  129 | 0.1363 | 0.0484 | -0.239 | 0.057 | 0.006 |

## Output Tables

- `reports/tables/edge_bootstrap_sensitivity/edge_inclusion_correlation_results.csv`
- `reports/tables/edge_bootstrap_sensitivity/edge_inclusion_correlation_comparison.csv`
- `reports/tables/edge_bootstrap_sensitivity/edge_nonedge_range_comparison.csv`
- `reports/tables/edge_bootstrap_sensitivity/equal_n_non_edge_resampling_results.csv`
- `reports/tables/edge_bootstrap_sensitivity/edge_nonedge_environment_comparison.csv`
- `reports/tables/edge_bootstrap_sensitivity/environment_correlations.csv`
- `reports/tables/edge_bootstrap_sensitivity/environment_adjusted_diversity_models.csv`
- `reports/tables/edge_bootstrap_sensitivity/bootstrap_variation_correlation_sensitivity.csv`
- `reports/tables/edge_bootstrap_sensitivity/analysis_dataset_with_bootstrap_fields.csv`

## Figures

- `reports/figures/edge_bootstrap_sensitivity/10m_spec_sa_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_mean_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_med_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_sd_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_rao_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_alpha_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_convex_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_hull3d_v_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_spec_spca_hull3d_a_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_sa_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_mean_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_med_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_sd_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_rao_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_alpha_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_convex_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_hull3d_v_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_spec_spca_hull3d_a_edge_highlight_scatter.png`
- `reports/figures/edge_bootstrap_sensitivity/10m_sa_entropy_bootstrap_cv_sensitivity.png`
- `reports/figures/edge_bootstrap_sensitivity/20m_sa_entropy_bootstrap_cv_sensitivity.png`
- `reports/figures/edge_bootstrap_sensitivity/equal_n_non_edge_resampling_primary_measures.png`
- `reports/figures/edge_bootstrap_sensitivity/edge_nonedge_environment_boxplots.png`

## Interpretation Notes

- A positive `delta r` means the correlation became stronger or more positive after edge quadrats were removed.
- A negative `delta r` means the correlation became weaker or more negative after edge quadrats were removed.
- `edge-only minus non-edge r` compares the edge-quadrat correlation directly against the non-edge correlation.
- In the equal-n resampling table, small tail probabilities mean the edge-only correlation is unusual relative to randomly selected non-edge subsets of the same size.
- Environmental adjusted models are screening diagnostics, not final spatial inference; the TRI variables are related topographic summaries and should be interpreted cautiously when fitted together.
- Bootstrap sensitivity applies directly to SA entropy only; standardized PCA metrics do not have bootstrap SD/CV fields in the current diagnostics table.
