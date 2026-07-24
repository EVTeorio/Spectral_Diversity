# Multiscale Spectral-Biodiversity Analysis

Last updated: 2026-07-10

## Current Direction

The active analysis has been revamped around the first research question: what is the relationship between spectral variation and phylogenetic/species diversity?

## Inputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- `RESEARCH_OBJECTIVES.md` for scientific framing.

## Primary Spectral Variation Measures

- `spec_spca_mean`: standardized PCA mean Euclidean distance.
- `spec_sa`: spectral angle entropy.

## Diversity Predictors

- Phylogenetic: `phy_faith`, `phy_rao`, `phy_afaith`.
- Species: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.

## Methods

- Each SV measure is independently paired with each diversity measure within each spatial scale.
- Each pair is summarized with Pearson `r`, `R2`, simple-regression `F`, F-test p-value, slope, intercept, and Spearman `r`.
- `sv_diversity_analysis_dataset.csv` includes `sa_all_pixels_sampled` and `spca_euclidean_all_pixels_sampled`, which flag whether each primary SV metric used all retained pixels for each quadrat.
- No environmental covariates or multivariable candidate-model ranking are included in this first relationship layer.
- Documented 10 m and 20 m edge quadrats are excluded from primary analysis; 50 m uses all quadrats.

## Outputs

- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`
- `reports/figures/multiscale_spectral_biodiversity/01_sv_diversity_pairwise_correlation_heatmap.png`
- `reports/figures/multiscale_spectral_biodiversity/02_sa_entropy_diversity_scatterplots.png`
- `reports/figures/multiscale_spectral_biodiversity/03_standardized_pca_distance_diversity_scatterplots.png`

## Result Size

- Pairwise correlation rows: 42
- Analysis dataset rows: 2580

## Superseded Material

The previous environment-adjusted AIC model-ranking workflow is superseded for the current stage. It can be revisited later as a second layer after the direct SV-diversity relationships are interpreted.
