# Task Report: Species And Phylogenetic Diversity Correlation Analysis

Date: 2026-08-17

## Objective

Compare species-diversity and phylogenetic-diversity metrics across 10 m, 20 m, and 50 m quadrats, and evaluate whether phylogenetic diversity adds explanatory value beyond species diversity for observed spectral variation.

## Outputs Created

- `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`
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
- `reports/figures/species_phylogenetic_correlation/`

## Result Size

- Species-phylogenetic pairwise correlation rows: 36
- All diversity metric pairwise correlation rows: 63
- Elevation-adjusted all-metric model rows: 63
- Incremental spectral-variation model rows: 72
- Species-composition cluster rows: 2580
- Named composition type rows: 15
- Species/individual moderation model rows: 126
- Significant moderation rows after BH adjustment: 52
- Residual divergence rows: 30960

## Notes

- Biodiversity-only comparisons use all quadrats with available species and phylogenetic diversity values.
- Spectral increment checks use all complete spectral-diversity cases and do not remove edge quadrats.
- Species-composition clusters are based on species presence/absence, not crown-overlap abundance.
- Composition type names are generated from the species codes most characteristic of each presence-based cluster.
- The `09_...composition_cluster...` scatterplot legends report panel-specific mean silhouette width for each composition type, calculated from the two metrics plotted in that panel.
- The `10_...individual_ramp...` scatterplot figures color quadrats by summed crown-overlap proportions, interpreted as crown-equivalent individuals rather than raw stem counts.
- The `11_...species_count_ramp...` scatterplot figures color quadrats by the number of species present, matching the presence-based richness definition.
- Species-count and crown-equivalent-individual moderation tests use standardized `y_metric ~ x_metric * moderator` models; significant result summaries exclude cases where species count is already the plotted richness metric.
- Species-composition clusters are descriptive, not formal community types.
