# Species And Phylogenetic Diversity Correlation Analysis

Date: 2026-08-17

## Purpose

This analysis compares four species-diversity measures with three phylogenetic-diversity measures across 10 m, 20 m, and 50 m quadrats. It asks how closely the phylogenetic metrics track species diversity, where they diverge, and whether adding phylogenetic diversity contributes additional explanatory variation in the observed spectral-heterogeneity data.

## Inputs

- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`
- `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv`
- `Quad_Values/Diversity_SHPs/plant_diversity_20m.csv`
- `Quad_Values/Diversity_SHPs/plant_diversity_50m.csv`

## Measures Compared

- Species-diversity measures: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.
- Phylogenetic-diversity measures: `phy_faith`, `phy_rao`, `phy_afaith`.
- Observed spectral-variation responses for incremental checks: `spec_spca_mean` and `spec_sa`.

## Methods

- Pairwise species-phylogenetic concordance was calculated within each scale using all quadrats with available species and phylogenetic diversity values.
- Each pair was summarized with Pearson correlation, linear-model R2, slope, intercept, F statistic, and Spearman correlation.
- A separate all-metric scatterplot set compares every unique pair among the four species-diversity metrics and three phylogenetic-diversity metrics within each quadrat scale.
- Edge-highlight scatterplot panels were created for 10 m and 20 m quadrats, with edge quadrats in red, non-edge quadrats in blue, and the regression line plus displayed statistics calculated from non-edge quadrats.
- Elevation-gradient scatterplot panels were created for each scale, with points colored from red at low mean elevation to green at high mean elevation. Each panel reports the incremental R2 and p-value for adding `env_elev` after the x-axis diversity metric.
- Individual-ramp scatterplot panels were created for each scale, with points colored by crown-equivalent individuals calculated as the sum of species crown-overlap proportions in each quadrat.
- Species-count scatterplot panels were created for each scale, with points colored by the number of species present in each quadrat.
- Incremental spectral-variation models compared `SV ~ species_metric`, `SV ~ phylo_metric`, and `SV ~ species_metric + phylo_metric` using all complete spectral-diversity cases; the reported incremental R2 is the additional variance explained by the phylogenetic metric after the species metric.
- To evaluate convergence or divergence around broadly similar species composition, quadrats were clustered within each scale from species presence/absence, not abundance. Phylogenetic residuals from `phylo_metric ~ species_metric` were summarized within those composition clusters.
- Mean silhouette width was calculated for each composition type within every all-metric scatterplot panel using the plotted two-metric space; larger positive values indicate stronger separation from the other composition types in that panel.
- Species count and crown-equivalent individuals were then tested as moderators of each all-metric relationship using standardized models of the form `y_metric ~ x_metric * moderator`. Significant results below use BH-adjusted interaction p-values and exclude species-count tests where species richness is already one of the plotted axes.

## Closest Species-Diversity Match For Each Phylogenetic Metric

| Scale | Phylogenetic metric | Closest species metric | n | r | R2 | p-value |
| --- | --- | --- | --- | --- | --- | --- |
| 10m | Faith's PD | Species richness | 2000 | 0.710 | 0.504 | 2.13e-306 |
| 10m | Phylogenetic Rao's Q | Simpson diversity | 2000 | 0.681 | 0.464 | 5.43e-273 |
| 10m | Abundance-weighted Faith's PD | Species richness | 2000 | 0.452 | 0.205 | 1.62e-101 |
| 20m | Faith's PD | Species richness |  500 | 0.727 | 0.528 | 3.23e-83 |
| 20m | Phylogenetic Rao's Q | Simpson diversity |  500 | 0.488 | 0.238 | 3.17e-31 |
| 20m | Abundance-weighted Faith's PD | Species richness |  500 | 0.236 | 0.056 | 8.78e-08 |
| 50m | Faith's PD | Species richness |   80 | 0.761 | 0.580 | 2.45e-16 |
| 50m | Phylogenetic Rao's Q | Species richness |   80 | 0.386 | 0.149 | 4.01e-04 |
| 50m | Abundance-weighted Faith's PD | Species evenness |   80 | -0.223 | 0.050 | 0.047 |

## Strongest Incremental Phylogenetic Contribution To Observed Spectral Variation

| Scale | Observed SV metric | Species metric already in model | Added phylogenetic metric | n | Species-only R2 | Combined R2 | Phylo incremental R2 | Partial p-value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | Std PCA mean distance | Species richness | Faith's PD | 1744 | 0.003 | 0.047 | 0.044 | 8.80e-19 |
| 10m | SA entropy | Simpson diversity | Phylogenetic Rao's Q | 1903 | 0.001 | 0.032 | 0.031 | 7.54e-15 |
| 20m | Std PCA mean distance | Species richness | Faith's PD |  436 | 0.001 | 0.141 | 0.140 | 7.04e-16 |
| 20m | SA entropy | Simpson diversity | Phylogenetic Rao's Q |  485 | 0.002 | 0.068 | 0.066 | 9.41e-09 |
| 50m | Std PCA mean distance | Simpson diversity | Phylogenetic Rao's Q |   74 | 0.000 | 0.223 | 0.223 | 2.43e-05 |
| 50m | SA entropy | Simpson diversity | Phylogenetic Rao's Q |   80 | 0.000 | 0.133 | 0.133 | 9.73e-04 |

## Largest Composition-Cluster Divergence Cases

| Scale | Cluster | Cluster size | Species metric | Phylogenetic metric | n | Mean abs residual z |
| --- | --- | --- | --- | --- | --- | --- |
| 20m | 20m_C5 | 90 | Species richness | Abundance-weighted Faith's PD | 90 | 1.104 |
| 50m | 50m_C1 | 26 | Shannon diversity | Phylogenetic Rao's Q | 26 | 1.091 |
| 20m | 20m_C5 | 90 | Shannon diversity | Abundance-weighted Faith's PD | 90 | 1.090 |
| 20m | 20m_C5 | 90 | Species evenness | Abundance-weighted Faith's PD | 90 | 1.090 |
| 50m | 50m_C1 | 26 | Species richness | Phylogenetic Rao's Q | 26 | 1.088 |
| 20m | 20m_C5 | 90 | Simpson diversity | Abundance-weighted Faith's PD | 90 | 1.081 |
| 20m | 20m_C3 | 68 | Shannon diversity | Faith's PD | 68 | 1.065 |
| 50m | 50m_C1 | 26 | Simpson diversity | Phylogenetic Rao's Q | 26 | 1.064 |
| 20m | 20m_C3 | 68 | Simpson diversity | Faith's PD | 68 | 1.059 |
| 20m | 20m_C3 | 68 | Species richness | Faith's PD | 68 | 1.056 |
| 50m | 50m_C1 | 26 | Species richness | Faith's PD | 26 | 1.048 |
| 20m | 20m_C3 | 68 | Species evenness | Faith's PD | 68 | 1.047 |

## Significant Species-Count Or Individual-Count Effects On Metric Relationships

| Scale | Relationship | Moderator | n | Interaction slope | Interaction delta R2 | Full delta R2 | Interaction q-value |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10m | Shannon diversity vs Simpson diversity | Number of species | 2000 | -0.107 | 0.017 | 0.037 | 4.99e-232 |
| 10m | Shannon diversity vs Species evenness | Number of species | 2000 | -0.153 | 0.036 | 0.283 | 9.37e-102 |
| 10m | Shannon diversity vs Simpson diversity | Crown-equivalent individuals | 2000 | -0.095 | 0.012 | 0.012 | 7.57e-93 |
| 10m | Phylogenetic Rao's Q vs Abundance-weighted Faith's PD | Crown-equivalent individuals | 2000 | 0.179 | 0.038 | 0.252 | 4.72e-82 |
| 10m | Shannon diversity vs Species evenness | Crown-equivalent individuals | 2000 | -0.186 | 0.045 | 0.048 | 2.71e-50 |
| 20m | Shannon diversity vs Simpson diversity | Number of species |  500 | -0.114 | 0.020 | 0.058 | 4.10e-43 |
| 10m | Species richness vs Species evenness | Crown-equivalent individuals | 2000 | -0.218 | 0.067 | 0.097 | 3.74e-36 |
| 10m | Species richness vs Simpson diversity | Crown-equivalent individuals | 2000 | -0.161 | 0.036 | 0.051 | 6.30e-36 |
| 20m | Simpson diversity vs Species evenness | Number of species |  500 | 0.189 | 0.041 | 0.195 | 4.74e-35 |
| 10m | Faith's PD vs Phylogenetic Rao's Q | Number of species | 2000 | -0.123 | 0.025 | 0.035 | 4.56e-34 |
| 10m | Simpson diversity vs Species evenness | Number of species | 2000 | -0.078 | 0.007 | 0.160 | 3.77e-31 |
| 10m | Faith's PD vs Abundance-weighted Faith's PD | Crown-equivalent individuals | 2000 | 0.108 | 0.015 | 0.215 | 7.32e-31 |

## Interpretation Notes

- Faith's PD is expected to converge most strongly with species richness because both are primarily presence-based.
- Phylogenetic Rao's Q and abundance-weighted Faith's PD can diverge from richness when common species differ in relatedness or when abundance is concentrated in phylogenetically similar or distinct lineages.
- A positive incremental R2 means the phylogenetic metric explains spectral-variation structure not captured by the paired species-diversity metric alone.
- Cluster residual summaries identify compositionally similar groups where phylogenetic values remain unusually high or low relative to the species-diversity expectation.
- A significant moderator interaction means the slope between two plotted diversity metrics changes as species count or crown-equivalent individuals increase.

## Figures

- `reports/figures/species_phylogenetic_correlation/01_species_phylogenetic_correlation_heatmap.png`
- `reports/figures/species_phylogenetic_correlation/02_phylogenetic_incremental_r2_heatmap.png`
- `reports/figures/species_phylogenetic_correlation/03_species_phylogenetic_scatter_grid.png`
- `reports/figures/species_phylogenetic_correlation/04_phylogenetic_divergence_by_metric_boxplot.png`
- `reports/figures/species_phylogenetic_correlation/05_all_diversity_metric_pairwise_scatter_10m.png`
- `reports/figures/species_phylogenetic_correlation/05_all_diversity_metric_pairwise_scatter_20m.png`
- `reports/figures/species_phylogenetic_correlation/05_all_diversity_metric_pairwise_scatter_50m.png`
- `reports/figures/species_phylogenetic_correlation/09_all_diversity_metric_pairwise_scatter_composition_cluster_10m.png`
- `reports/figures/species_phylogenetic_correlation/09_all_diversity_metric_pairwise_scatter_composition_cluster_20m.png`
- `reports/figures/species_phylogenetic_correlation/09_all_diversity_metric_pairwise_scatter_composition_cluster_50m.png`
- `reports/figures/species_phylogenetic_correlation/10_all_diversity_metric_pairwise_scatter_individual_ramp_10m.png`
- `reports/figures/species_phylogenetic_correlation/10_all_diversity_metric_pairwise_scatter_individual_ramp_20m.png`
- `reports/figures/species_phylogenetic_correlation/10_all_diversity_metric_pairwise_scatter_individual_ramp_50m.png`
- `reports/figures/species_phylogenetic_correlation/11_all_diversity_metric_pairwise_scatter_species_count_ramp_10m.png`
- `reports/figures/species_phylogenetic_correlation/11_all_diversity_metric_pairwise_scatter_species_count_ramp_20m.png`
- `reports/figures/species_phylogenetic_correlation/11_all_diversity_metric_pairwise_scatter_species_count_ramp_50m.png`
- `reports/figures/species_phylogenetic_correlation/06_all_diversity_metric_pairwise_scatter_edge_highlight_10m.png`
- `reports/figures/species_phylogenetic_correlation/06_all_diversity_metric_pairwise_scatter_edge_highlight_20m.png`
- `reports/figures/species_phylogenetic_correlation/07_all_diversity_metric_pairwise_scatter_elevation_gradient_10m.png`
- `reports/figures/species_phylogenetic_correlation/07_all_diversity_metric_pairwise_scatter_elevation_gradient_20m.png`
- `reports/figures/species_phylogenetic_correlation/07_all_diversity_metric_pairwise_scatter_elevation_gradient_50m.png`
- `reports/figures/species_phylogenetic_correlation/08_species_presence_composition_type_heatmap_10m.png`
- `reports/figures/species_phylogenetic_correlation/08_species_presence_composition_type_heatmap_20m.png`
- `reports/figures/species_phylogenetic_correlation/08_species_presence_composition_type_heatmap_50m.png`

## Output Tables

- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_pairwise_correlations.csv`
- `reports/tables/species_phylogenetic_correlation/all_diversity_metric_pairwise_correlations.csv`
- `reports/tables/species_phylogenetic_correlation/elevation_adjusted_all_diversity_metric_models.csv`
- `reports/tables/species_phylogenetic_correlation/phylogenetic_incremental_sv_models.csv`
- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_metric_summary.csv`
- `reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv`
- `reports/tables/species_phylogenetic_correlation/species_composition_type_names.csv`
- `reports/tables/species_phylogenetic_correlation/species_composition_type_presence_frequencies.csv`
- `reports/tables/species_phylogenetic_correlation/species_composition_scatterplot_silhouette.csv`
- `reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv`
- `reports/tables/species_phylogenetic_correlation/species_individual_moderated_metric_relationships.csv`
- `reports/tables/species_phylogenetic_correlation/significant_species_individual_moderated_metric_relationships.csv`
- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_residual_divergence.csv`
- `reports/tables/species_phylogenetic_correlation/species_composition_cluster_divergence_summary_presence_based.csv`
