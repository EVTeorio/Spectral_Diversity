# Validation: Quadrat 800_a Alpha-Hull Figure

## Validation Steps

- Ran the figure script with R 4.2.3 using the standard user library path.
- Recalculated retained illuminated pixels, unique PC1-PC2 points, convex-hull area, and alpha-hull area from `Quad_Spectra/10m_smooth_5nm/800_a`.
- Compared recalculated values with `Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`.
- Visually inspected the exported PNG.

## Results

| Check | Recalculated | Existing table | Status |
|---|---:|---:|---|
| Retained illuminated pixels | 3,976 | 3,976 | Pass |
| Unique PC1-PC2 points | 3,900 | 3,900 | Pass |
| Convex-hull area | 6,679.321799829 | 6,679.321799826 | Pass |
| Alpha-hull area | 628.063050156 | 628.063050156 | Pass |

## Conclusion

The exported figure represents the actual vector-normalized standardized PCA alpha-hull calculation for quadrat `800_a` with `alpha = 1`. The recalculated alpha-hull and convex-hull areas match the existing spectral heterogeneity metric table to numerical precision.
