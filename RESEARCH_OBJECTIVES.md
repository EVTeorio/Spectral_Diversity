# Research Objective and Current Analysis Specification

## Project

Evaluate whether canopy spectral heterogeneity from drone hyperspectral imagery tracks species and phylogenetic diversity in the Paint Rock Forest Dynamics Plot across quadrat scales.

## Current Analytical Focus

The active analysis uses both original and vector-normalized global PCA spaces to describe spectral heterogeneity. The current priority metric families are:

- mean Euclidean distance of pixels from the quadrat centroid in PCA space;
- equal-pixel-weight Rao's Q with squared Euclidean distance in PC1-PC3 space;
- 2D alpha-hull area in PCA PC1-PC2 space.

Both PCA bases and their corresponding distance and alpha-hull metrics are retained. Vector normalization is an important brightness-reduced comparison, not a replacement for the original PCA basis.

Current biodiversity responses are species richness, Shannon diversity, Simpson diversity, evenness, Faith's phylogenetic diversity, phylogenetic Rao's Q, and abundance-weighted Faith's phylogenetic diversity. These are retained in the plant-diversity outputs and combined quadrat tables.

Primary questions are whether these relationships vary by scale and whether phylogenetic diversity tracks spectral heterogeneity differently from species diversity.

## Inputs and Scales

Current confirmed base spectra are stored in `Quad_Spectra/10m/`, `20m/`, `25m/`, and `50m/`. The current resampled inputs are the corresponding `*_smooth_5nm/` folders. The 10 m folder was restored from a data repository on October 15; the dated verification compares its vector-normalized PC1-PC3 mean distance and Rao's Q against saved values. Missing `*_smooth/` intermediates were deliberately removed because they are reproducible.

The current established scales are 10 m, 20 m, and 50 m. The active expansion is 25 m; all-pixel PC1-PC3 mean distance, Rao's Q, and PC1-PC2 alpha = 0.5 areas are complete for its 312 prepared rasters in both saved PCA bases. The user handles exclusions manually. The authorized alpha = 0.5 passes at 10 m, 20 m, and 25 m are finished, with only original-PCA 10 m 1000_c unresolved and deferred for later inspection. A user-authorized manual arc-orientation repair resolved 20 m original-PCA 413, with separate provenance. The user also reviewed and authorized orientation repairs for 10 m 1402_b, 1019_a, and 1816_c; their corrected values retain explicit manual provenance. Processing stopped before 50 m. 33 m and 100 m are deferred.

## Canonical Current Data Products

- Spectral summaries and both PCA basis objects: `Quad_Values/Spectral_diversitySHPs/`
- 25 m distance/Rao-only summary: `Quad_Values/Spectral_diversitySHPs/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv` (312 prepared rasters, 320 grid rows; no exclusions applied)
- Current alpha = 0.5 areas: `Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_*m_all_pixels.csv` for 10 m, 20 m, and 25 m (both bases; explicit full-grid statuses)
- Species and phylogenetic diversity summaries: `Quad_Values/Diversity_SHPs/`
- Environmental covariates: `Quad_Values/Enviro_SHPs/`
- Analysis-ready joined tables: `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`

The original-PCA fields in the joined tables are `spec_pca_mean`, `spec_rao_q`, and `spec_alpha`; the vector-normalized equivalents are `spec_spca_mean`, `spec_spca_rao`, and `spec_spca_alpha`. Future joined tables should preserve these names or clearly map equivalent fields.

## Archive and Provenance

Spectral-angle entropy, convex hulls, 3D hull metrics, bootstrap uncertainty products, and sample-size sensitivity outputs remain useful for provenance and possible later review, but are not primary forward metrics. They should be preserved compactly, not treated as active results by default.

## Interpretation Boundaries

Analyses remain observational. Spatial autocorrelation and environmental context require sensitivity treatment before causal or final inferential language. The current manuscript materials in `Documents/` are protected and are outside routine repository consolidation.
