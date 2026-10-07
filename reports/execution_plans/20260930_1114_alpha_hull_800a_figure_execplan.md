# ExecPlan: 10 m Quadrat 800_a Alpha-Hull Figure

## Objective

Create a manuscript-ready image using 10 m quadrat `800_a` to show the difference between the convex hull and the alpha hull in vector-normalized PCA space, using the same alpha-hull calculation applied in the spectral heterogeneity workflow.

## Requested Task

Generate an image that depicts the actual alpha-hull area calculated for quadrat `800_a`, alongside the convex hull for comparison.

## Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `reports/directory_map.md`
- `reports/data_dictionary.md`
- `tools/generate_alpha_vs_convex_hull_10m_example.R`
- `tools/export_10m_vector_pca_scores_for_alpha_examples.R`
- `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`

## Relevant Files

- `Quad_Spectra/10m_smooth_5nm/800_a`
- `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds`
- `Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv`

## Proposed Changes

- Add a targeted R plotting script under `scripts/visuals/`.
- Export a PNG comparison figure under `Documents/Tables and Figures/`.
- Export a small metadata table under `reports/tables/figure_sources/`.

## Expected Modifications

- `scripts/visuals/create_800a_alpha_vs_convex_hull_true_figure.R`
- `Documents/Tables and Figures/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated.png`
- `reports/tables/figure_sources/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated_metadata.csv`

## Validation Plan

- Use R 4.2.3 and the standard R 4.2 user library.
- Recalculate retained pixel count, unique PC1-PC2 point count, convex-hull area, and alpha-hull area.
- Compare recalculated values against the existing 10 m spectral heterogeneity summary row for `800_a`.
- Confirm the PNG file exists and has nonzero file size.

## Risks

- Alpha-hull plotting from `alphahull` produces arcs rather than a simple polygon, so the displayed alpha hull will emphasize the true calculated boundary and reported area rather than attempting an approximate filled polygon that could misrepresent discontinuities.
- If package dependencies are unavailable in the standard R library, figure generation would require dependency installation approval.
