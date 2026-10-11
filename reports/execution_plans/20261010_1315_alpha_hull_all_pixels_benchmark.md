# Objective

Benchmark all-retained-pixel 2D alpha-hull area calculations for a small representative subset of quadrats across 10 m, 20 m, and 50 m scales, using R `alphahull` and multiple alpha parameter values.

# Requested Task

Run an R-only benchmark before revisiting Python. Use 8 selected 10 m quadrats, 4 selected 20 m quadrats, and 2 selected 50 m quadrats. Compare resulting alpha-hull values with abundance-weighted Faith's PD.

# Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/directory_map.md`
- `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`
- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`

# Relevant Files

- `Quad_Spectra/10m_smooth_5nm/`
- `Quad_Spectra/20m_smooth_5nm/`
- `Quad_Spectra/50m_smooth_5nm/`
- `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`
- `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds`

# Proposed Changes

- Add a focused R benchmark script under `scripts/3_Analysis/`.
- Select representative non-excluded quadrats by retained-pixel count within scale.
- Recompute PC1-PC2 alpha-hull areas from all retained pixels for original and vector-normalized PCA bases.
- Evaluate multiple alpha values, including the current production value `alpha = 1`.
- Write benchmark tables and an analysis report under `reports/`.

# Expected Modifications

- New script: `scripts/3_Analysis/alpha_hull_all_pixels_benchmark.R`
- New report: `reports/analysis/20261010_alpha_hull_all_pixels_benchmark.md`
- New tables under `reports/tables/alpha_hull_all_pixels_benchmark/`
- New task and validation notes under `reports/tasks/` and `reports/validation/`
- Updates to `reports/project_state.md`, `reports/data_dictionary.md`, and `reports/directory_map.md`

# Validation Plan

- Confirm required R packages are available.
- Confirm selected quadrats have matching spectra, PCA values, and abundance-weighted Faith's PD.
- Confirm output tables are non-empty and contain all requested scales.
- Confirm correlation table reports Pearson `r` versus `phy_afaith`.
- Confirm failures/timeouts are explicitly represented in output fields.

# Risks

- Full all-pixel `alphahull` may be slow or fail for 50 m quadrats with very large retained-pixel counts.
- The benchmark subset is small, so correlation estimates are diagnostic rather than final inference.
- Alpha-hull area may be sensitive to alpha value, point density, duplicated points, and fragmented components.
