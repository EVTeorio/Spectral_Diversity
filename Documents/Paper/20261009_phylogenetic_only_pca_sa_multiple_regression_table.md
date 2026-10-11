# Phylogenetic-Only Multiple Regression Models for PCA and Spectral Angle Heterogeneity

OLS model form for each scale and spectral metric:

`spectral_metric ~ Faith's PD + phylogenetic Rao's Q + abundance-weighted Faith's PD`

| Spectral_metric | 10m n | 10m R2 | 10m adj R2 | 10m p | 20m n | 20m R2 | 20m adj R2 | 20m p | 50m n | 50m R2 | 50m adj R2 | 50m p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Spectral angle entropy | 1903 | 0.014 | 0.013 | 5.88e-06 | 485 | 0.044 | 0.038 | 8.28e-05 | 80 | 0.120 | 0.086 | 0.020 |
| Original PCA alpha hull | 1741 | 0.005 | 0.004 | 0.025 | 436 | 0.048 | 0.041 | 9.63e-05 | 74 | 0.177 | 0.141 | 0.003 |
| Original PCA mean distance | 1744 | 0.004 | 0.002 | 0.074 | 436 | 0.027 | 0.020 | 0.008 | 74 | 0.145 | 0.108 | 0.012 |
| Original PCA Rao's Q | 1744 | 0.001 | -0.001 | 0.568 | 436 | 0.017 | 0.010 | 0.065 | 74 | 0.133 | 0.096 | 0.018 |
| Vector-normalized PCA alpha hull | 1744 | 0.038 | 0.037 | 1.05e-14 | 436 | 0.152 | 0.146 | 2.25e-15 | 74 | 0.289 | 0.259 | 2.44e-05 |
| Vector-normalized PCA mean distance | 1744 | 0.020 | 0.019 | 9.13e-08 | 436 | 0.096 | 0.090 | 1.61e-09 | 74 | 0.222 | 0.189 | 4.99e-04 |
| Vector-normalized PCA Rao's Q | 1744 | 0.002 | 0.000 | 0.279 | 436 | 0.021 | 0.014 | 0.028 | 74 | 0.126 | 0.088 | 0.024 |

Source: `reports/tables/pca_sa_phylogenetic_only_multiple_regression/pca_sa_phylogenetic_only_full_model_summary.csv`.

Notes: spectral angle entropy is included as the angular non-PCA comparison metric. Original PCA refers to the brightness-sensitive PCA basis; vector-normalized PCA refers to the PCA basis built after vector-normalizing spectra. Models use complete cases from the current combined quadrat tables. Significance is from the full-model F test.
