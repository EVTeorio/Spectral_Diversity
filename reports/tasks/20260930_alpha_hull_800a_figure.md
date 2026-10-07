# Task Report: Quadrat 800_a Alpha-Hull Figure

## Objective

Create a figure using the actual 10 m quadrat `800_a` spectral PCA point cloud to compare convex-hull area with the alpha-hull area used in the current spectral heterogeneity metrics.

## Actions Performed

- Reviewed `CODEX_AGENT_GUIDELINES.md` and switched R dependency handling to the standard R 4.2 user library.
- Confirmed `alphahull` is available from `C:/Users/PaintRock/AppData/Local/R/win-library/4.2`.
- Reviewed the active spectral heterogeneity workflow and confirmed alpha-hull area is calculated with `alphahull::ahull(..., alpha = 1)` and `alphahull::areaahull()`.
- Added `scripts/visuals/create_800a_alpha_vs_convex_hull_true_figure.R`.
- Generated a two-panel figure comparing the convex hull and true alpha-hull boundary for `800_a` in vector-normalized standardized PCA PC1-PC2 space.
- Saved a metadata CSV comparing recalculated values with the existing metric table values.

## Outputs

- `Documents/Tables and Figures/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated.png`
- `reports/tables/figure_sources/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated_metadata.csv`

## Notes

The figure preserves the alpha-hull object's discontinuous zero-radius components, drawn smaller than the default `alphahull` plotting method so the main occupied-space boundary remains legible.
