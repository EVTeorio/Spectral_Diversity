# Moran's I for Spectral Metrics, Biodiversity Metrics, and Spectral-Biodiversity Residuals

Moran's I values for the primary vector-normalized PCA spectral metrics, biodiversity metrics, and residuals from priority spectral-biodiversity relationships at 10, 20, and 50 m quadrat grains.

## Spectral Metrics

| Metric | 10 m Moran's I | 20 m Moran's I | 50 m Moran's I |
|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | 0.417 | 0.318 | 0.221 |
| Vector-normalized PCA mean Euclidean distance | 0.281 | 0.263 | 0.214 |
| Vector-normalized PCA spectral Rao's Q | 0.178 | 0.174 | 0.112 |

## Biodiversity Metrics

| Metric | 10 m Moran's I | 20 m Moran's I | 50 m Moran's I |
|---|---:|---:|---:|
| Phylogenetic Rao's Q | 0.295 | 0.368 | 0.373 |
| Abundance-weighted Faith's PD | 0.349 | 0.367 | 0.283 |
| Shannon diversity | 0.284 | 0.333 | 0.214 |

## Spectral-Biodiversity Relationship Residuals

| Spectral metric | Biodiversity metric | 10 m residual Moran's I | 20 m residual Moran's I | 50 m residual Moran's I |
|---|---|---:|---:|---:|
| Vector-normalized PCA alpha-hull area | Phylogenetic Rao's Q | 0.394 | 0.228 | 0.051 |
| Vector-normalized PCA alpha-hull area | Abundance-weighted Faith's PD | 0.385 | 0.202 | 0.068 |
| Vector-normalized PCA alpha-hull area | Shannon diversity | 0.416 | 0.318 | 0.219 |
| Vector-normalized PCA mean Euclidean distance | Phylogenetic Rao's Q | 0.262 | 0.192 | 0.073 |
| Vector-normalized PCA mean Euclidean distance | Abundance-weighted Faith's PD | 0.262 | 0.179 | 0.101 |
| Vector-normalized PCA mean Euclidean distance | Shannon diversity | 0.281 | 0.263 | 0.212 |
| Vector-normalized PCA spectral Rao's Q | Phylogenetic Rao's Q | 0.175 | 0.151 | 0.041 |
| Vector-normalized PCA spectral Rao's Q | Abundance-weighted Faith's PD | 0.176 | 0.151 | 0.080 |
| Vector-normalized PCA spectral Rao's Q | Shannon diversity | 0.177 | 0.175 | 0.104 |

Source: `reports/tables/final_research_direction/spatial_moran_diagnostics.csv`.

Spatial weights used eight nearest neighbors. Moran's I permutation tests used 199 permutations in the source workflow.
