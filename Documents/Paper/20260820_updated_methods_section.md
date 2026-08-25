# Updated Methods Section

Created: 2026-08-20

## Methods

### Study Design

This study evaluated the Spectral Variation Hypothesis across three quadrat grains in the Paint Rock Forest Dynamics Plot: 10 m, 20 m, and 50 m. The workflow paired UAV hyperspectral canopy spectra with field-derived biodiversity metrics calculated for spatially matched quadrats. Analyses were conducted within each scale so that spectral heterogeneity, species diversity, phylogenetic diversity, transformation behavior, and spatial structure could be compared across quadrat grain.

The final PCA-centered workflow used three forward spectral heterogeneity metrics as the main PCA set: vector-normalized standardized PCA alpha-hull area, vector-normalized standardized PCA mean Euclidean distance, and vector-normalized standardized PCA Rao's Q. Raw PCA analogs were retained as brightness-sensitive comparisons. Convex hull and 3D hull metrics were calculated and archived as diagnostic metric-screening outputs, but they were not carried forward as primary manuscript metrics.

### Field Tree And Quadrat Data

Field tree records were read from `PR_tree_DL.csv`. Taxonomic information was read from `51sp_taxanomy.csv`. Quadrat geometries were read from `Quad_Scale_SHPs/PR_10m.shp`, `Quad_Scale_SHPs/PR_20m.shp`, and `Quad_Scale_SHPs/PR_50m.shp`.

Tree records were filtered to retain individuals with `DBH.2024 >= 200`, crown position classes 3, 4, or 5, cluster status `A` or `R`, nonmissing current UTM coordinates, nonmissing positive crown width, and species codes represented in the taxonomy table. The `COOB2` taxon was excluded from the taxonomy table before phylogeny construction. Crown radii were calculated as one half of `cw_m_2025`. Tree points were projected to EPSG:26916 and buffered by crown radius to create circular crown polygons.

For each scale, tree crown polygons were intersected with quadrat polygons. The area of each crown-quadrat intersection was divided by the total crown area to estimate the proportion of each crown overlapping each quadrat. Crown-overlap proportions were summed by species within each quadrat to create a canopy-weighted species abundance matrix. These summed crown-overlap values represent crown-equivalent abundance rather than raw stem count.

### Biodiversity Metrics

Species and phylogenetic diversity metrics were calculated from the crown-overlap abundance matrix for each quadrat.

Species richness, `sp_rich`, was calculated as the number of species with positive crown-overlap abundance. Shannon diversity, `sp_shannon`, was calculated from proportional species abundances as `-sum(p * log(p))`. Simpson diversity, `sp_simpson`, was calculated as `1 - sum(p^2)`. Species evenness, `sp_even`, was calculated as Shannon diversity divided by `log(richness)` when richness was greater than one; quadrats with richness less than or equal to one were assigned evenness of zero.

A phylogeny was generated from genus, species, and family information using `V.PhyloMaker2` with `GBOTB.extended.TPL`, `nodes.info.1.TPL`, and scenario 3. Faith's phylogenetic diversity, `phy_faith`, was calculated as the summed branch length of the subtree connecting species present in a quadrat. Phylogenetic Rao's Q, `phy_rao`, was calculated from proportional crown-overlap abundance and pairwise cophenetic distances among species in the quadrat. Abundance-weighted Faith's PD, `phy_afaith`, was calculated by multiplying branch lengths by abundance weights assigned to terminal branches from the corresponding species crown-overlap abundances and summing the weighted branch lengths.

The all-scale biodiversity workflow was implemented in `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`. Outputs were written to `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv`, `Quad_Values/Diversity_SHPs/plant_diversity_20m.csv`, and `Quad_Values/Diversity_SHPs/plant_diversity_50m.csv`, with matching shapefiles.

### Hyperspectral Imagery And Spectral Preprocessing

Hyperspectral canopy data were organized as quadrat-level raster stacks for 10 m, 20 m, and 50 m quadrats. The active spectral inputs were the smoothed and 5 nm resampled rasters in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.

Initial quadrat spectra were smoothed using the project Savitzky-Golay smoothing workflow and then resampled to 5 nm spectral spacing. Sidecar files with extensions such as `.hdr`, `.aux`, `.xml`, `.enp`, and `.sta` were excluded from raster-file iteration. Spectra were cleaned by removing rows with missing values, zero-only spectra, nonfinite values, and spectra not retained by the illumination mask.

### Shadow Masking And Illumination Screening

Shadow and sunlit training polygons were stored in `Shadow_vs_sunlit_SHP/Shadow_vs_sunlit_Rev.shp`. The shadow-screening diagnostic extracted pixel values from the training polygons across the 20 m smoothed 5 nm raster set and evaluated wavelength-specific discrimination using ROC curves and AUC. For each candidate wavelength, the polygon label was treated as the class response and band reflectance was treated as the predictor. The optimal threshold was identified from ROC coordinates using the best-threshold criterion.

The active spectral workflow retained sunlit pixels using the 563 nm band and a reflectance threshold of `0.0305476`. Pixels with 563 nm reflectance greater than this threshold were retained. This mask was applied before spectral angle entropy, raw PCA metric, and vector-normalized PCA metric calculations. The shadow/sunlit AUC diagnostic was implemented in `scripts/1_Data Processing/Sunlit_vs_Shadow_AUC.R`; the retained mask was applied in the spectral heterogeneity scripts through the same 563 nm threshold rule.

### PCA Basis Construction

Raw and vector-normalized PCA bases were constructed from current 10 m smoothed 5 nm rasters after shadow masking. The PCA basis construction sampled up to 450 retained illuminated pixels from each 10 m raster, using deterministic seed logic based on quadrat identifiers. The raw PCA basis was fit to the retained spectra without vector normalization. The vector-normalized standardized PCA basis was fit after dividing each retained spectrum by its Euclidean norm.

Both PCA bases used 121 spectral bands and retained three principal components for quadrat-level metric calculations. PCA objects and variance-explained tables were stored in `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`, `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`, `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds`, and `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv`.

The active PCA exclusion policy was `ten_meter_footprint_pca_450_pixels_20260707`. Manual atmospheric/cloud exclusions were applied to documented 20 m and 50 m quadrats; matching 10 m exclusions were derived from the excluded 20 m quadrats. Excluded quadrats were retained as rows with missing PCA metric values and method fields indicating manual exclusion.

### Spectral Angle Entropy

Spectral angle entropy, `spec_sa`, was calculated from sunlit, shadow-masked, smoothed 5 nm spectra. Pairwise spectral angles were calculated from retained spectra and summarized with Shannon entropy over binned angle values. For small enough quadrats, all retained pixel pairs were used. For larger quadrats, entropy was calculated as the mean of 70 bootstrap iterations using up to 5,000 retained pixels per iteration.

Spectral angle entropy outputs were generated by `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R`. Scale-specific summary files were written to `Quad_Values/10m_SA_entropy_smooth_masked_5nm_summary.csv`, `Quad_Values/20m_SA_entropy_smooth_masked_5nm_summary.csv`, and `Quad_Values/50m_SA_entropy_smooth_masked_5nm_summary.csv`. Bootstrap long and wide tables were also written for each scale and retained as quality-control outputs.

### PCA-Based Spectral Heterogeneity Metrics

PCA-based spectral heterogeneity metrics were calculated from retained pixels projected into raw PCA and vector-normalized standardized PCA space. Metrics were calculated for each quadrat and each scale using `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`.

For mean PCA distance, retained pixel scores in PC1-PC3 space were centered by the quadrat centroid in PCA space. The metric was calculated as the mean Euclidean distance from each retained pixel to that centroid. Median distance and standard deviation of distances were also calculated and archived as supporting summaries.

Spectral Rao's Q was calculated in PC1-PC3 space using equal pixel weights and squared Euclidean distance. The implemented value used the equivalent centroid formula `2 * mean(squared_radius)`, where `squared_radius` is the squared distance from each retained pixel to the quadrat centroid in PCA space.

PCA alpha-hull area was calculated from retained pixel scores in PC1-PC2 space using the `alphahull` package with alpha set to 1. When the retained point count exceeded configured limits, points were sampled according to the alpha-hull workflow settings before hull calculation. The workflow archived method fields indicating whether all pixels, sampled pixels, fallback sampling, or an exclusion condition was used.

Raw PCA metrics were written to `Quad_Values/*_PCA_metrics_smooth_masked_5nm_summary.csv`. Vector-normalized standardized PCA metrics were written to `Quad_Values/*_standardized_PCA_metrics_smooth_masked_5nm_summary.csv`. Combined spectral heterogeneity tables were written to `Quad_Values/*_spectral_heterogeneity_smooth_masked_5nm_summary.csv` and `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_*_smooth_masked_5nm_summary.csv`, with matching shapefiles in `Quad_Values/Spectral_diversitySHPs/`.

### Environmental And Illumination Context Variables

Environmental summaries were calculated for each quadrat scale and joined into the combined analysis tables. Mean elevation, `env_elev`, came from the leaf-off digital terrain model.

Retained-pixel brightness variables were calculated from the same shadow-masked quadrat rasters used for spectral heterogeneity. The brightness summary included mean 563 nm reflectance and regional mean reflectance in blue, green, red, and near-infrared wavelength regions. These variables were stored in `reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv` and used in metric-driver analyses.

### Combined Quadrat Analysis Tables

Current spectral heterogeneity, species diversity, phylogenetic diversity, environmental summaries, and quadrat centroids were joined into one root-level table for each scale using `scripts/3_Analysis/combine_quadrat_analysis_tables.R`. The combined tables were `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.

Combined table centroids were taken from the plant-diversity shapefiles and stored as `center_x` and `center_y` in `NAD83 / UTM zone 16N`. Pixel-count, pair-count, bootstrap replicate, method, exclusion, geometry metadata, and per-species composition fields were excluded from the root-level combined tables. Shortened metric names were documented in `reports/combined_quadrat_variable_guide.md`.

### Primary Spectral-Biodiversity Correlation Analysis

All spectral-biodiversity relationships were evaluated within scale using complete cases for each metric pair. The all-metric scatter analysis used 6 spectral heterogeneity measures: raw PCA alpha hull, Rao's Q, mean distance; and standardized PCA alpha hull, Rao's Q, mean distance. Seven biodiversity metrics were evaluated: Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, species richness, Shannon diversity, Simpson diversity, and species evenness.

For each scale, spectral metric, and biodiversity metric, a simple linear model was fit with the spectral metric as the response and the biodiversity metric as the predictor. Pearson correlation, `R2`, F statistic, F-test p-value, slope, intercept, and complete-case sample size were saved. Each plot panel used the corresponding complete-case data for that metric pair. The all-metric workflow did not apply edge filtering; upstream missingness and documented PCA exclusions were preserved.

This analysis was implemented in `scripts/3_Analysis/spectral_biodiversity_all_metric_scatter_analysis.R`. Tables were written to `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_dataset.csv` and `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_relationships.csv`. Figures were written to `reports/figures/spectral_biodiversity_all_metrics/`, with one scale-column figure per spectral metric.

### Biodiversity Metric Concordance Analysis

Relationships among biodiversity metrics were evaluated across all available biodiversity quadrats using `scripts/3_Analysis/species_phylogenetic_correlation_analysis.R`. The analysis compared four species-diversity measures (`sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`) with three phylogenetic-diversity measures (`phy_faith`, `phy_rao`, `phy_afaith`) within each scale.

The workflow calculated pairwise correlations, all-diversity metric summaries, species presence/absence composition clusters, species-composition type names, species-composition silhouette summaries, crown-equivalent individual totals, and species-count totals. Scatterplot panels were generated with points colored by species-composition type, crown-equivalent individuals, or present species count. Standardized moderation models of the form `y_metric ~ x_metric * moderator` were also fit for species count and crown-equivalent individuals; species-count moderation tests were excluded when species richness was already one of the plotted axes. Benjamini-Hochberg adjustment was applied to interaction p-values.

Outputs were written to `reports/tables/species_phylogenetic_correlation/`, `reports/figures/species_phylogenetic_correlation/`, and `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`.

### Spectral Metric Concordance And Driver Analysis

Relationships among spectral heterogeneity metrics were evaluated using `scripts/3_Analysis/spectral_heterogeneity_relationship_analysis.R`. The workflow used the scale-specific combined quadrat tables and complete cases for each spectral metric pair. The analysis produced pairwise scatterplots and correlation tables for spectral heterogeneity metrics across 10 m, 20 m, and 50 m.

Additional plotting layers colored the same spectral metric relationships by mean elevation, species presence/absence composition type, crown-equivalent individuals, present species count, mean retained-pixel 563 nm brightness, and retained-pixel blue, green, red, and near-infrared brightness. Spectral metric driver relationships were summarized using pairwise complete-case correlations between focal spectral metrics and environmental or illumination-context variables.

Outputs were written to `reports/tables/spectral_heterogeneity_relationships/`, `reports/figures/spectral_heterogeneity_relationships/`, and `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`.

### Raw PCA Versus Vector-Normalized PCA Diagnostics

Two diagnostics were used to document how the raw PCA and vector-normalized standardized PCA bases differed.

First, 50 m mean-reflectance diagnostics sampled 250 valid sunlit pixels from each 50 m quadrat for both PCA bases. For each sampled pixel, mean reflectance across spectral bands was calculated and compared with PC1, PC2, and PC3 scores. Pearson correlations were calculated at the pixel level and after averaging to one mean point per quadrat. This workflow was implemented in `scripts/3_Analysis/pc_mean_reflectance_correlation_50m_both.R`, with outputs in `reports/tables/pc1_mean_reflectance_correlation/` and reports in `reports/analysis/20260706_50m_pc1_mean_reflectance_correlation.md` and `reports/analysis/20260707_50m_pc_mean_reflectance_correlation.md`.

Second, PCA loading diagnostics compared PC1 and PC2 loading structure between raw and vector-normalized standardized PCA. PC1 loadings were oriented so the mean loading was positive and compared with a uniform brightness vector of `1 / sqrt(n_bands)`. Wavelengths with the largest absolute departures from the brightness vector were identified for PC1. Wavelengths with the largest absolute PC2 loadings were identified for PC2. Contiguous 5 nm wavelengths meeting each rule were collapsed into spectral regions. This analysis was implemented in `scripts/3_Analysis/pca_loading_spectral_region_analysis.R`, with outputs in `reports/tables/pca_loading_spectral_regions/`, `reports/figures/pca_loading_spectral_regions/`, and `reports/analysis/20260707_pca_loading_spectral_regions.md`.

### Transformation Analyses

A broad transformation screen was run using `scripts/3_Analysis/final_research_direction_analysis.R`. Candidate transformations included identity, log, `log1p`, square root, inverse, `power_0.25`, `power_0.5`, `power_2`, Box-Cox, and Yeo-Johnson transformations. Transformations were applied to spectral and biodiversity axes where mathematically valid. For each transformed metric pair, the workflow saved Pearson correlation, Spearman correlation, `R2`, linear-model statistics, Shapiro-Wilk residual p-value when applicable, and correlation between absolute residuals and fitted values.

The final reported transformation sensitivity was restricted to the biodiversity-squared axis comparison. A second transformation workflow, `scripts/3_Analysis/spectral_biodiversity_all_metric_scatter_analysis.R`, created two power-2 scenarios for every spectral-biodiversity metric pair: biodiversity squared with original spectral values, and biodiversity squared with spectral values also squared. The retained manuscript text reports the biodiversity-squared/original-spectral scenario. Transformation relationship tables were written to `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_power2_transform_relationships.csv`. Corresponding figures were written to `reports/figures/spectral_biodiversity_all_metrics/` and stacked phylogenetic transformation figures were written to `reports/figures/spectral_biodiversity_phylogenetic_transform_stacks/`.

### Spatial Autocorrelation Diagnostics

Spatial diagnostics were conducted for priority variables and model residuals using `scripts/3_Analysis/final_research_direction_analysis.R`. Priority spectral metrics were standardized PCA alpha hull, standardized PCA mean distance, and standardized PCA Rao's Q. Priority biodiversity metrics were phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity.

Moran's I was calculated within scale using quadrat centroid coordinates and eight-nearest-neighbor spatial weights. Permutation tests used 199 permutations. Variable-level Moran's I was calculated for the priority spectral and biodiversity metrics. Residual Moran's I was calculated after fitting scale-specific linear models of each priority spectral metric as a function of each priority biodiversity metric. A semivariogram summary was also calculated for variables and residuals.

Spatial diagnostic tables were written to `reports/tables/final_research_direction/spatial_moran_diagnostics.csv` and `reports/tables/final_research_direction/spatial_variogram_summary.csv`. The spatial diagnostic figure was written to `reports/figures/final_research_direction/02_spatial_moran_i_priority_variables_and_residuals.png`.

### Metric Driver And Context Analyses

Metric-driver analyses were generated as part of `scripts/3_Analysis/final_research_direction_analysis.R` and supporting spectral/biodiversity relationship scripts. Spectral metrics were compared with elevation, TRI variables, 563 nm brightness, regional brightness summaries, retained-pixel brightness variables, species-composition type, present species count, and crown-equivalent individuals where available. Biodiversity metrics were compared with present species count, crown-equivalent individuals, species-composition clusters, and other biodiversity metrics.

Driver relationship tables were written to `reports/tables/final_research_direction/metric_driver_relationships.csv` and related context-specific tables in `reports/tables/spectral_heterogeneity_relationships/` and `reports/tables/species_phylogenetic_correlation/`.

### Edge, Bootstrap, And Sample-Size Sensitivity Outputs

Edge-quadrat and bootstrap sensitivity analyses were retained as supporting workflow outputs. Edge analyses compared relationships with and without documented 10 m and 20 m edge quadrats; 50 m quadrats were retained because no separate 50 m edge rule was documented. Equal-sample-size non-edge resampling and environmental screening models were generated as supporting checks. These analyses were implemented in `scripts/3_Analysis/edge_and_bootstrap_sensitivity_sv_diversity.R`, with outputs under `reports/tables/edge_bootstrap_sensitivity/`, `reports/figures/edge_bootstrap_sensitivity/`, and `reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md`.

Bootstrap variation diagnostics for spectral angle entropy were generated using `scripts/3_Analysis/bootstrap_variation_analysis.R` and `scripts/3_Analysis/sa_bootstrap_context_relationship_analysis.R`. These workflows summarized bootstrap standard deviation, coefficient of variation, standard error, and confidence-interval fields for spectral angle entropy, then related these uncertainty fields to environmental, brightness, and biodiversity-context variables.

Sample-size sensitivity analyses were run separately from the final main PCA values. Spectral angle entropy sample-size sensitivity was implemented in `scripts/3_Analysis/sa_entropy_sample_size_effects.R`. PCA mean distance, spectral Rao's Q, and alpha-hull area sample-size sensitivity were implemented in `scripts/3_Analysis/spectral_metric_sample_size_effects.R` for both raw and vector-normalized PCA bases. Fixed 4,000-pixel outputs from these workflows were treated as sensitivity outputs rather than as final main PCA metric values.

### Figure And Table Generation

The updated analysis outputs were generated from the current combined quadrat tables and supporting diagnostic tables. The main output families were:

- All spectral-biodiversity scatterplots and relationship tables in `reports/figures/spectral_biodiversity_all_metrics/` and `reports/tables/spectral_biodiversity_all_metrics/`.
- Phylogenetic transformation stack figures in `reports/figures/spectral_biodiversity_phylogenetic_transform_stacks/`.
- Standardized PCA volume/hull biodiversity comparison figure in `reports/figures/spectral_biodiversity_all_metrics/`.
- Spectral metric relationship figures and brightness/elevation context figures in `reports/figures/spectral_heterogeneity_relationships/`.
- Biodiversity metric relationship, composition-cluster, species-count, and crown-equivalent individual figures in `reports/figures/species_phylogenetic_correlation/`.
- Transformation, metric-driver, and spatial-diagnostic tables in `reports/tables/final_research_direction/`.
- PCA mean-reflectance and loading diagnostic outputs in `reports/tables/pc1_mean_reflectance_correlation/` and `reports/tables/pca_loading_spectral_regions/`.

All correlation and regression summaries used complete cases for the relevant metric pair and scale. Missingness from upstream atmospheric exclusions, shadow-masked spectra, insufficient retained pixels, or hull computation conditions was preserved rather than imputed.
