# Combined Quadrat Analysis Tables

Last updated: 2026-06-24

## Task

Combined current spectral heterogeneity, species and phylogenetic diversity, environmental/topographic values, and quadrat center coordinates into one root-level CSV for each quadrat scale.

## Outputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- `combined_quadrat_variable_guide.docx`
- `reports/combined_quadrat_variable_guide.md`
- `reports/validation/20260624_combined_quadrat_analysis_tables_validation.md`

## Notes

- Pixel-count, pair-count, bootstrap replicate, method, exclusion, and geometry metadata fields were excluded.
- Per-species composition columns were excluded; original composition values remain available in `Quad_Values/Diversity_SHPs/plant_diversity_*m.csv` if needed later.
- Species diversity summaries use the `sp_` prefix.
- Spectral columns use the `spec_` prefix.
- Environmental/topographic columns use the `env_` prefix.
- Spectral gaps remain `NA` where current raster-derived spectral summaries were not available or where PCA-dependent values were excluded upstream.

## Validation Summary

| Scale | Rows | Columns | Duplicate quad IDs | Missing center X | Missing center Y | Missing spectral SA | Missing spectral PCA mean | Missing elevation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10m | 2000 | 24 | 0 | 0 | 0 | 97 | 256 | 0 |
| 20m |  500 | 24 | 0 | 0 | 0 | 15 |  64 | 0 |
| 50m |   80 | 24 | 0 | 0 | 0 |  0 |   6 | 0 |
