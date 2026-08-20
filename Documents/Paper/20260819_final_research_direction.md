# Final Research Direction

Created: 2026-08-19

## Working Manuscript Focus

This paper should be framed as a multiscale test of whether UAV hyperspectral spectral heterogeneity captures phylogenetic diversity in a temperate forest dynamics plot better than conventional species-diversity summaries.

The final manuscript should emphasize a cautious, correlation-based claim:

> Drone hyperspectral spectral heterogeneity shows modest but consistent positive relationships with phylogenetic diversity, especially at broader spatial grains, while standard species-diversity metrics show weak, near-zero, or inconsistent relationships.

The strongest current direction is not that spectral heterogeneity is a universal biodiversity proxy. The defensible claim is that spectral heterogeneity appears more aligned with phylogenetic structure than with species diversity alone, and that this relationship is scale dependent.

## Primary Research Question

Can UAV hyperspectral spectral heterogeneity serve as a spatial proxy for field-derived phylogenetic diversity in the Paint Rock Forest Dynamics Plot, and does that relationship strengthen from 10 m to 20 m to 50 m quadrat grains?

## Final Objectives

1. Quantify spectral heterogeneity from sunlit, shadow-masked, smoothed 5 nm UAV hyperspectral spectra across 10 m, 20 m, and 50 m quadrats.
2. Compare spectral heterogeneity against field-derived species and phylogenetic diversity metrics calculated from crown-overlap-weighted tree census data.
3. Test whether phylogenetic diversity metrics explain spectral heterogeneity better than conventional species-diversity metrics.
4. Use species-phylogenetic concordance results to explain why phylogenetic metrics can add signal beyond richness, Shannon diversity, Simpson diversity, and evenness.
5. Present edge, bootstrap, elevation, brightness, and metric-stability analyses as sensitivity or supporting layers rather than as the central manuscript claim.

## Recommended Final Thesis

The paper should argue that hyperspectral spectral heterogeneity is most informative when biodiversity is represented as phylogenetic composition rather than as species diversity alone. The relationship is weak at 10 m, clearer at 20 m, and strongest at 50 m, suggesting that broader grains better integrate canopy spectral mixtures, crown overlap, species composition, and evolutionary relatedness.

The final interpretation should remain careful because the 50 m results have fewer quadrats, the analyses are mostly pairwise correlations, and spatial autocorrelation has not yet been incorporated as a final inferential layer.

## Primary Spectral Metrics

Use these as the main spectral heterogeneity metrics:

- `spec_spca_mean`: standardized PCA mean Euclidean distance. This is the cleanest primary PCA metric because vector normalization reduces broad brightness dominance and it has already been central in manuscript drafts.
- `spec_spca_alpha`: standardized PCA alpha-hull area. The 2026-08-19 all-metric analysis shows this is currently the strongest spectral-biodiversity metric, especially at 50 m.
- `spec_sa`: spectral angle entropy. This provides an independent angular spectral-heterogeneity metric and has the strongest bootstrap uncertainty documentation.

Treat these as supporting or supplemental spectral metrics:

- `spec_spca_rao`: standardized PCA spectral Rao's Q.
- `spec_pca_mean`, `spec_alpha`, and `spec_rao_q`: raw PCA counterparts, useful as brightness-sensitive comparisons.
- Convex hull and 3D hull summaries: supplemental diagnostics only.

## Primary Biodiversity Metrics

Use these as the main biodiversity metrics:

- `phy_rao`: phylogenetic Rao's Q. This is the most consistent phylogenetic response, especially for standardized PCA mean distance and SA entropy.
- `phy_afaith`: abundance-weighted Faith's PD. This pairs strongly with standardized PCA alpha hull and standardized PCA mean distance, especially at 20 m and 50 m.
- `sp_shannon`: main species-diversity contrast, because prior paper notes already identify it as the clearest conventional diversity comparison.

Use these as supporting metrics:

- `phy_faith`: important for explaining presence-based phylogenetic diversity and its strong convergence with richness.
- `sp_rich`, `sp_simpson`, and `sp_even`: include in all-metric figures/tables or supplement to show that the phylogenetic signal is not simply a generic species-diversity effect.

## Core Result Story

The results should be organized around four linked findings.

1. Spectral-biodiversity relationships are phylogenetically stronger than species-diversity relationships.

The July 10 primary analysis found the strongest pairings for both primary spectral measures were phylogenetic at every scale. Species richness, Shannon diversity, Simpson diversity, and evenness were weak, near zero, or negative in the primary non-edge analysis.

2. Relationships strengthen with spatial grain.

The clearest primary standardized PCA result is at 50 m: `spec_spca_mean` versus `phy_rao` has r = 0.444 and R2 = 0.197. The 20 m signal is moderate, led by `spec_spca_mean` versus `phy_afaith` with r = 0.291 and R2 = 0.085 in the non-edge primary analysis. The 10 m signal is positive but weak.

3. The newest all-metric analysis identifies standardized PCA alpha hull as a lead metric.

The August 19 all-metric scatter analysis found the strongest absolute relationships at 50 m, led by `spec_spca_alpha` with `phy_afaith` (r = 0.497, R2 = 0.247) and `phy_rao` (r = 0.487, R2 = 0.238). At 20 m, `spec_spca_alpha` also performed strongly with `phy_afaith` (r = 0.384, R2 = 0.147).

4. Phylogenetic metrics add information beyond species metrics.

The species-phylogenetic concordance analysis shows Faith's PD closely tracks richness, but phylogenetic Rao's Q and abundance-weighted Faith's PD diverge more strongly from conventional species metrics. Incremental models show phylogenetic metrics add explanatory signal for spectral variation after paired species metrics, with the largest gains at 50 m for phylogenetic Rao's Q.

## Recommended Main Figures

1. Study system and workflow figure.
2. Primary all-metric spectral-biodiversity figure for `spec_spca_alpha`.
3. Primary all-metric spectral-biodiversity figure for `spec_spca_mean`.
4. Spectral angle entropy biodiversity figure as independent confirmation.
5. Species-phylogenetic concordance heatmap or scatter-grid figure.
6. Sensitivity figure showing edge and bootstrap patterns, likely as supplement unless space allows.

## Recommended Main Tables

1. Quadrat coverage and complete-case counts by scale.
2. Metric definitions for spectral heterogeneity, species diversity, and phylogenetic diversity.
3. Top spectral-biodiversity relationships by scale and metric.
4. Species-phylogenetic concordance and incremental phylogenetic contribution summary.

Place full pairwise correlation tables, edge comparisons, bootstrap-CV sensitivity, brightness/elevation diagnostics, and all supplemental spectral metric relationships in supplementary material.

## Sensitivity Results Placement

Use the sensitivity analyses to support caution rather than to redirect the paper.

- Edge quadrats: edge-only correlations can be stronger, but edge and non-edge range differences complicate interpretation. Main figures should either use primary non-edge filtering for 10 m and 20 m or clearly flag edge quadrats.
- Bootstrap uncertainty: SA entropy low-CV subsets retain or strengthen phylogenetic correlations; high-CV subsets weaken. This supports using SA bootstrap fields as quality-control context.
- Elevation: elevation is negatively associated with spectral variation, especially at 20 m. Environmental-adjusted models suggest phylogenetic diversity retains some incremental signal, but this should remain a screening result unless spatial/environmental models are finalized.
- Brightness: raw PCA metrics remain useful comparisons, but standardized PCA metrics should lead because the standardized PCA workflow reduces brightness dominance.

## What To Avoid Chasing Before This Paper Draft

- Do not center the manuscript on every spectral metric. Lead with standardized PCA alpha hull, standardized PCA mean distance, and spectral angle entropy.
- Do not make environmental model ranking the main result. Keep environmental adjustment as sensitivity unless a final spatial model is completed.
- Do not overbuild the species-composition clustering narrative. Use it to explain divergence between species and phylogenetic metrics, not as a separate community-type paper.
- Do not treat 50 m p-values as final inference without noting the small sample size and need for spatial diagnostics.

## Final Output Checklist

- Final manuscript outline updated to the current lead metrics.
- Results text for `spec_spca_alpha`, `spec_spca_mean`, and `spec_sa`.
- Figure inventory with final main-versus-supplement assignments.
- Table inventory with top results and complete-case counts.
- Short sensitivity paragraph for edge, bootstrap uncertainty, elevation, and brightness.
- Spatial autocorrelation or spatial sensitivity check before final inferential language.
- MDPI Remote Sensing-compatible abstract, keywords, conclusions, data availability statement, and cover-letter draft.

## Source Markdown Reviewed

This direction is based on the current paper workspace notes, project-state documentation, data guides, analysis reports, task reports, and validation notes, especially:

- `Documents/Paper/20260730_working_outline.md`
- `Documents/Paper/20260730_paper_workspace_index.md`
- `Documents/Paper/20260803_writing_style_notes.md`
- `Documents/Paper/20260803_markdown_context_digest.md`
- `Documents/Paper/20260803_remote_sensing_mdpi_author_requirements.md`
- `RESEARCH_OBJECTIVES.md`
- `reports/project_state.md`
- `reports/data_dictionary.md`
- `reports/directory_map.md`
- `reports/combined_quadrat_variable_guide.md`
- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`
- `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`
- `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`
- `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`
- `reports/analysis/20260817_bootstrap_metric_stat_relationship_analysis.md`
- `reports/analysis/20260817_sa_bootstrap_context_relationship_analysis.md`
- `reports/analysis/20260819_spectral_biodiversity_all_metric_scatter_analysis.md`
- Matching 2026-08-17 and 2026-08-19 task and validation notes.
