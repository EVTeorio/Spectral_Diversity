# Data Dictionary

Last updated: 2026-07-06

This dictionary inventories current datasets visible in the repository. Large raster and hyperspectral products were documented from file metadata and project context; row/column counts are listed where lightweight tabular inspection was possible. Unnamed first columns in CSV outputs are treated as likely index columns unless later script review shows otherwise.

Path note: as of 2026-07-06, the live derived-output directory visible in the working tree is `Quad_Values/`. Producing scripts and documentation now use this path for derived spectral, diversity, environmental, topographic, and index outputs.

## Tabular Datasets

### `PR_tree_DL.csv`

- Dataset Name: Paint Rock tree census and trait table
- Description: Tree-level census, spatial, taxonomic, DBH, crown, biomass, and trait fields used to derive biodiversity metrics.
- Source Location: `PR_tree_DL.csv`
- Row Count: 28,774
- Column Count: 69
- Primary Keys: `OBJECTID`, `tag`, `StemTag`, `uniqueID` are candidate identifiers; confirm before joins.
- Date Fields: `date.measured`, `ExactDate`, `CreationDate`, `EditDate`
- Numeric Fields: `OBJECTID`, coordinates, quadrat, tag fields, DBH fields, crown metrics, biomass fields, wood specific gravity, growth fields
- Categorical Fields: species codes, status fields, notes, family/genus/species, common name, taxonomic level, height class, native/clonal flags, references
- Missing Value Summary: Missingness is substantial in older census, notes, crown, lean, and some trait/growth fields; completely missing columns include `ave.dbh`, `max.dbh`, `ave.grownth`, and `max.growth`.
- Last Updated: 2026-03-09

| Variable Group | Type | Description |
|---|---|---|
| `OBJECTID`, `tag`, `StemTag`, `uniqueID` | identifier | Tree and stem identifiers |
| `new.gx`, `new.gy`, `UTMX_CURRENT`, `UTMY_CURRENT`, `quadrat` | numeric/spatial | Tree coordinates and quadrat assignment |
| `sp`, `Family`, `Genus`, `Species`, `Common.Name` | categorical | Species and taxonomy fields |
| `previous_dbh`, `DBH.2024`, `DBH.2024_in` | numeric | Diameter measurements |
| `crown.position`, `height.class`, `cw_*` | mixed | Crown position and crown width fields |
| `BA2020`, `BA2025`, `agb2020`, `agb2025` | numeric | Basal area and aboveground biomass estimates |
| `status.2024`, `previous_status`, `census_status`, `codes.2024` | categorical | Census status and annotation fields |

### `51sp_taxanomy.csv`

- Dataset Name: 51 species taxonomy table
- Description: Species code lookup with species, genus, family, order, and clade.
- Source Location: `51sp_taxanomy.csv`
- Row Count: 51
- Column Count: 6
- Primary Keys: `sp_code`
- Date Fields: none
- Numeric Fields: none detected
- Categorical Fields: `sp_code`, `species`, `genus`, `family`, `order`, `clade`
- Missing Value Summary: none detected
- Last Updated: 2026-03-09

| Variable | Type | Description |
|---|---|---|
| `sp_code` | categorical | Species code |
| `species` | categorical | Scientific species name |
| `genus` | categorical | Genus |
| `family` | categorical | Family |
| `order` | categorical | Order |
| `clade` | categorical | Higher-level phylogenetic grouping |

### `species_summary.csv`

- Dataset Name: Species abundance summary
- Description: Species-level count summary, likely derived from tree census filtering.
- Source Location: `species_summary.csv`
- Row Count: 51
- Column Count: 2
- Primary Keys: `sp`
- Date Fields: none
- Numeric Fields: `n_individuals`
- Categorical Fields: `sp`
- Missing Value Summary: none detected
- Last Updated: 2026-03-23

| Variable | Type | Description |
|---|---|---|
| `sp` | categorical | Species code |
| `n_individuals` | numeric | Count of individuals for each species |

### `spec_shannaon_DIV_with Geometry.csv`

- Dataset Name: Spectral and Shannon diversity with geometry
- Description: Quadrat-level table combining Shannon diversity and spectral metric with geometry.
- Source Location: `spec_shannaon_DIV_with Geometry.csv`
- Row Count: 500
- Column Count: 5
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: unnamed index column, `Name`, `shannon`, `spctrl_`
- Categorical Fields: `geometry`
- Missing Value Summary: none detected
- Last Updated: 2026-02-05

| Variable | Type | Description |
|---|---|---|
| unnamed first column | numeric | Likely row/index column |
| `Name` | numeric | Quadrat identifier |
| `shannon` | numeric | Shannon diversity metric |
| `spctrl_` | numeric | Spectral metric field, exact derivation should be confirmed |
| `geometry` | geometry text | Geometry representation |

### `Quad_Values/20m_spectral_sp.csv`

- Dataset Name: 20 m spectral, biodiversity, species, and environmental table
- Description: Main integrated 20 m analysis table combining spectral entropy/PCA metrics, phylogenetic and species-diversity metrics, species abundance columns, geometry, and environmental variables.
- Source Location: `Quad_Values/20m_spectral_sp.csv`
- Row Count: 500
- Column Count: 69
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: spectral metrics, diversity metrics, species-code abundance columns, dominance/living-tree/elevation fields
- Categorical Fields: `geometry`
- Missing Value Summary: spectral metrics have 64 missing values each: `SA_entropy_smooth_masked_711`, `SA_entropy_smooth_masked`, `SA_entropy_smooth`, `global_PCA_20m_masked_5nm`, `SA_entropy_20m_masked_5nm`, `SumBandEntropy_20m_masked_5nm`
- Last Updated: 2026-04-20

| Variable Group | Type | Description |
|---|---|---|
| `Name` | numeric identifier | 20 m quadrat identifier |
| `SA_entropy_*`, `SumBandEntropy_*` | numeric | Spectral entropy metrics |
| `global_PCA_20m_masked_5nm` | numeric | Global PCA spectral heterogeneity metric |
| `AfaithPD`, `faithPD`, `raoPD` | numeric | Phylogenetic diversity metrics |
| `shnnn_d`, `smpsn_d`, `rchnss_`, `evnnss_` | numeric | Species diversity metrics |
| species-code columns such as `CAGL8`, `FAGR`, `LITU` | numeric | Species abundance or weighted species matrix fields |
| `dmnnt_v`, `avg_lvt`, `elvtn_r` | numeric | Dominance/living-tree/elevation related covariates |
| `geometry` | geometry text | Quadrat geometry |

### `Quad_Values/20m_SA_smooth_masked.csv`

- Dataset Name: 20 m masked smoothed spectral angle entropy
- Description: Quadrat-level spectral entropy table using smoothed and masked spectra.
- Source Location: `Quad_Values/20m_SA_smooth_masked.csv`
- Row Count: 436
- Column Count: 3
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: unnamed index column, `Name`, `spectral_entropy`
- Categorical Fields: none detected
- Missing Value Summary: none detected
- Last Updated: 2026-04-01

| Variable | Type | Description |
|---|---|---|
| unnamed first column | numeric | Likely row/index column |
| `Name` | numeric | Quadrat identifier |
| `spectral_entropy` | numeric | Spectral angle entropy |

### `Quad_Values/20m_SA_smooth_masked_7_11.csv`

- Dataset Name: 20 m masked smoothed spectral angle entropy, 7/11 variant
- Description: Quadrat-level spectral entropy table using a smoothed/masked parameter variant.
- Source Location: `Quad_Values/20m_SA_smooth_masked_7_11.csv`
- Row Count: 436
- Column Count: 3
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: unnamed index column, `Name`, `spectral_entropy`
- Categorical Fields: none detected
- Missing Value Summary: none detected
- Last Updated: 2026-04-01

| Variable | Type | Description |
|---|---|---|
| unnamed first column | numeric | Likely row/index column |
| `Name` | numeric | Quadrat identifier |
| `spectral_entropy` | numeric | Spectral angle entropy |

### `Quad_Values/20m_SA_entrop_boot_results.csv`

- Dataset Name: 20 m spectral angle entropy bootstrap results, 10 iterations
- Description: Quadrat-level bootstrap values for spectral angle entropy.
- Source Location: `Quad_Values/20m_SA_entrop_boot_results.csv`
- Row Count: 436
- Column Count: 18
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: `Name`, `n_pixels`, `boot_1` through `boot_10`, and bootstrap summary fields
- Categorical Fields: none detected
- Missing Value Summary: none detected
- Last Updated: 2026-05-06

| Variable Group | Type | Description |
|---|---|---|
| `Name` | numeric identifier | Quadrat identifier |
| `n_pixels` | numeric | Number of pixels sampled |
| `boot_1` to `boot_10` | numeric | Bootstrap replicate entropy values |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary statistics |

### `Quad_Values/20m_SA_entrop_boot100_results.csv`

- Dataset Name: 20 m spectral angle entropy bootstrap results, 100 iterations
- Description: Current expanded bootstrap values for spectral angle entropy.
- Source Location: `Quad_Values/20m_SA_entrop_boot100_results.csv`
- Row Count: 436
- Column Count: 108
- Primary Keys: `Name`
- Date Fields: none detected
- Numeric Fields: `Name`, `n_pixels`, `boot_1` through `boot_100`, and bootstrap summary fields
- Categorical Fields: none detected
- Missing Value Summary: none detected
- Last Updated: 2026-05-07

| Variable Group | Type | Description |
|---|---|---|
| `Name` | numeric identifier | Quadrat identifier |
| `n_pixels` | numeric | Number of pixels sampled |
| `boot_1` to `boot_100` | numeric | Bootstrap replicate entropy values |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary statistics |

### Current spectral angle entropy outputs from `_smooth_5nm`

- Dataset Name: Current per-scale spectral angle entropy from smoothed 5 nm spectra
- Description: Per-quadrat spectral heterogeneity tables generated by `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R` from current smoothed 5 nm spectra. For most quadrats, `spectral_entropy` is the mean of 70 bootstrap iterations using up to 5,000 retained pixels per iteration; large bootstrap iterations use 10,000 sampled pixel pairs. The small subset with pair counts below the configured exact threshold uses all retained pixel pairs. Therefore, these outputs should not be described as exhaustive all-pair calculations except for rows with `method = "exact_all_pixels"`.
- Source Location: `Quad_Values/*_SA_entropy_smooth_masked_5nm_*.csv`
- Row Count: summary CSVs contain 1,909 rows for 10 m, 485 rows for 20 m, and 80 rows for 50 m; long bootstrap CSVs contain 133,630 rows for 10 m, 33,950 rows for 20 m, and 5,600 rows for 50 m
- Column Count: summary CSVs contain 12 columns; wide bootstrap CSVs contain 72 columns, including 70 bootstrap replicate columns
- Primary Keys: `quad_id`
- Date Fields: none detected
- Numeric Fields: `n_pixels`, `n_pairs`, `spectral_entropy`, `exact_entropy`, `boot_mean`, `boot_sd`, `boot_cv`, bootstrap replicate fields
- Categorical Fields: `quad_id`, `method`
- Missing Value Summary: 10 m summary has 6 missing `spectral_entropy` values from insufficient valid pixels; 20 m and 50 m summaries have no missing `spectral_entropy` values. Output shapefiles have 97 missing `spec_ent` values for 10 m, 15 for 20 m, and 0 for 50 m.
- Last Updated: 2026-06-18

| Variable | Type | Description |
|---|---|---|
| `quad_id` | identifier | Quadrat identifier matched to the appropriate scale: 10 m `sub_id`, 20 m `Name`, or 50 m `Name` |
| `n_pixels` | numeric | Count of sunlit, valid spectra after shadow masking and NA/zero filtering |
| `n_pairs` | numeric | Number of possible pairwise spectral angles among retained pixels; for `bootstrap_mean` rows this is metadata, not the number of pairs actually evaluated |
| `method` | categorical | `exact_all_pixels`, `bootstrap_mean`, or `insufficient_pixels` |
| `spectral_entropy` | numeric | Primary spectral heterogeneity value; all-retained-pixel entropy for exact rows, otherwise the mean of 70 bootstrap iterations |
| `exact_entropy` | numeric | All-pixel spectral angle entropy when pair counts are below the configured exact threshold |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary fields used when all-pixel entropy is not computationally reasonable |
| `boot_1` ... `boot_n` | numeric | Wide-format bootstrap replicate values, where `n` is controlled by `SA_N_BOOT` |

### Current combined spectral heterogeneity outputs from `_smooth_5nm`

- Dataset Name: Current per-scale combined spectral heterogeneity outputs
- Description: Exclusion-aware spectral heterogeneity tables generated by `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`. These combine existing SA entropy with PCA-dependent metrics. The current PCA basis excludes manually documented atmospheric/cloud-affected quadrats before PCA sampling and leaves PCA-dependent values missing for those quadrats. Earlier PCA-dependent outputs from 2026-06-22 created before this exclusion policy should be disregarded.
- Source Location: `Quad_Values/*_spectral_heterogeneity_smooth_masked_5nm_summary.csv` and `Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_*_smooth_masked_5nm.*`
- Row Count: combined CSVs contain 1,909 rows for 10 m, 485 rows for 20 m, and 80 rows for 50 m
- Column Count: not re-counted after all joins; includes PCA-dependent metric fields and prefixed SA entropy fields
- Primary Keys: `quad_id`
- Date Fields: none detected
- Numeric Fields: `pca_euclidean_mean`, `pca_euclidean_median`, `pca_euclidean_sd`, `alpha_hull_area`, `pca_hull_volume_3d`, `pca_hull_area_3d`, `rao_q_pca`, `sa_entropy`, SA bootstrap fields
- Categorical Fields: `quad_id`, `scale`, `metric_method`, `manual_excluded`, `alpha_hull_method`, `pca_hull_3d_method`, `sa_method`
- Missing Value Summary: 10 m has 165 missing PCA Euclidean/Rao/3D hull values from manual exclusions and 168 missing alpha-hull values from those exclusions plus 3 alpha failures; 20 m has 49 missing PCA-dependent values from manual exclusions; 50 m has 6 missing PCA-dependent values from manual exclusions. SA entropy missingness is unchanged from the SA entropy workflow.
- Last Updated: 2026-06-22

| Variable | Type | Description |
|---|---|---|
| `quad_id`, `scale` | identifier | Quadrat identifier and analysis scale |
| `manual_excluded` | logical | TRUE when the quadrat is in the documented atmospheric/cloud exclusion set for PCA-dependent metrics |
| `metric_method` | categorical | `all_pixels` for included quadrats or `manual_excluded` for excluded quadrats |
| `pca_euclidean_mean` | numeric | Mean distance of all retained pixels from the quadrat centroid in global PC1-PC3 space |
| `alpha_hull_area` | numeric | Alpha-hull area in global PC1-PC2 space; may use deterministic sampling for large quadrats |
| `rao_q_pca` | numeric | Rao's Q using equal pixel weights and squared Euclidean distance in global PC1-PC3 space |
| `pca_hull_volume_3d`, `pca_hull_area_3d` | numeric | Supplemental convex hull volume and surface area in global PC1-PC3 space; not a true 3D alpha hull |
| `sa_*` fields | numeric/categorical | Prefixed fields imported from the current SA entropy summary CSVs; use `sa_method` to distinguish exact all-pixel entropy from bootstrap means |

### `reports/tables/bootstrap_variation/*.csv`

- Dataset Name: Bootstrap variation quality-control tables
- Description: Derived validation tables comparing within-quadrat bootstrap variation with between-quadrat spectral heterogeneity variation for 10 m, 20 m, and 50 m spectral entropy outputs.
- Source Location: `reports/tables/bootstrap_variation/`
- Row Count: 3 rows in `bootstrap_variation_scale_summary.csv`; 2,422 rows in `bootstrap_variation_quadrat_diagnostics.csv`; 45 rows in `bootstrap_variation_top_unstable_quadrats.csv`; 45 rows in `bootstrap_variation_widest_ci_quadrats.csv`
- Column Count: 42 columns in scale summary; 24 columns in quadrat diagnostics; 24 columns in top unstable quadrats; 24 columns in widest-CI quadrats
- Primary Keys: `scale` for scale summary; `scale` + `quad_id` for quadrat-level diagnostics
- Date Fields: none detected
- Numeric Fields: bootstrap SD/CV, between-scale variance summaries, standard errors, bootstrap-mean confidence intervals, pixel counts, pair counts, outlier replicate counts
- Categorical Fields: `scale`, `quad_id`, `method`, `stability_class`
- Missing Value Summary: exact all-pixel and insufficient-pixel quadrats are excluded from quadrat-level bootstrap diagnostics; no missing values expected in the scale summary
- Last Updated: 2026-06-18

| Variable Group | Type | Description |
|---|---|---|
| `scale`, `quad_id`, `method` | identifier/categorical | Analysis scale, quadrat identifier, and spectral entropy calculation method |
| `between_sd`, `between_var`, `between_iqr` | numeric | Between-quadrat spectral entropy variation at each scale |
| `mean_entropy`, `median_entropy`, `entropy_skewness` | numeric | Scale-level distribution summaries for final spectral heterogeneity values |
| `within_sd_*`, `within_cv_*` | numeric | Within-quadrat bootstrap variation summaries |
| `boot_se_*`, `boot_ci_*` | numeric | Standard error and 95 percent confidence interval summaries for the reported bootstrap mean |
| `within_to_between_var_ratio`, `icc_like_between_fraction` | numeric | Diagnostics comparing within-bootstrap variation to between-quadrat variation |
| `quads_cv_gt_005`, `quads_cv_gt_010` | numeric | Counts of quadrats exceeding 5% and 10% bootstrap CV thresholds |
| `outlier_reps` | numeric | Number of bootstrap replicates outside the quadrat-level 1.5 IQR fences |

### `reports/tables/sample_size_effects/sa_entropy/*.csv`

- Dataset Name: SA entropy sample-size effects experiment
- Description: Bootstrap sensitivity experiment generated by `scripts/3_Analysis/sa_entropy_sample_size_effects.R`. The workflow uses 32 selected 10 m quadrat spectra rasters, 16 selected 20 m rasters, and 8 selected 50 m rasters while retaining the original six pilot quadrats. It compares spectral angle entropy across scale-aware sample-size rules: 1%, 2%, and 3% of retained pixels capped at 5,000 pixels; fixed 1,250 and 4,000 pixels at all scales; fixed 2,000 and 3,000 pixels for 10 m and 20 m; and fixed 6,000 and 8,000 pixels for 50 m. Each quadrat x sample-rule combination uses 50 bootstrap iterations. When a rule resolves to 100% of retained pixels, the full-pixel entropy is calculated once and repeated across the 50 output rows so the full-pixel condition has zero artificial bootstrap variation. Outputs are stored under an `sa_entropy` spectral-variation-type folder so parallel experiments for other spectral metrics can be added cleanly.
- Source Location: `reports/tables/sample_size_effects/sa_entropy/`
- Row Count: 392 rows in `sa_entropy_sample_size_design.csv`; 19,600 rows in `sa_entropy_sample_size_boot_long.csv`; 392 rows in `sa_entropy_sample_size_summary.csv`
- Column Count: 18 columns in design; 21 columns in long bootstrap results; 29 columns in summary
- Primary Keys: `scale`, `quad_id`, `sample_rule`; plus `bootstrap_iter` for long bootstrap results
- Date Fields: none detected
- Numeric Fields: retained pixel counts, sample sizes, sample fractions, bootstrap entropy values, bootstrap means/SD/CV/SE, 95 percent confidence intervals, and differences from the fixed-4,000 rule
- Categorical Fields: `scale`, `quad_id`, `sample_rule`, `rule_type`, `rule_axis_label`, `sample_label`, `plot_sample_label`, `applicable_scales`, `pair_method`
- Missing Value Summary: no missing bootstrap entropy values expected; fixed 4,000 is capped to 3,976 pixels for 10 m quadrat `800_a`, which is the only current 100% retained-pixel condition in the expanded design
- Last Updated: 2026-07-04

| Variable Group | Type | Description |
|---|---|---|
| `scale`, `quad_id`, `sample_rule`, `rule_type`, `sample_label` | identifier/categorical | Quadrat scale and sample-size rule, with fixed-count labels showing the realized percent of each quadrat in parentheses |
| `n_pixels`, `sample_size`, `sample_fraction` | numeric | Retained pixel pool and number/fraction of pixels sampled per bootstrap iteration |
| `bootstrap_iter`, `spectral_entropy` | numeric | Bootstrap iteration index and resulting SA entropy value |
| `entropy_mean`, `entropy_sd`, `entropy_cv`, `entropy_se` | numeric | Mean and variability of the 50 bootstrap entropy values for each quadrat x sample-rule combination |
| `ci_low`, `ci_high`, `ci_half_width` | numeric | 95 percent confidence interval around the 50-iteration mean |
| `delta_from_fixed_4000` | numeric | Difference between each sample rule mean and the same quadrat's fixed-4,000-pixel mean |

### `reports/tables/sample_size_effects/{pca_mean_distance,spectral_rao_q,alpha_hull_area}/*.csv`

- Dataset Name: PCA-derived spectral metric sample-size effects experiments
- Description: Bootstrap sensitivity experiments generated by `scripts/3_Analysis/spectral_metric_sample_size_effects.R`. These reuse the finalized SA entropy quadrat/sample design so PCA mean distance, spectral Rao's Q, and alpha-hull area are evaluated on the same 32 selected 10 m quadrats, 16 selected 20 m quadrats, and 8 selected 50 m quadrats. Each quadrat x sample-rule combination uses 50 bootstrap iterations. Each replicate samples retained pixels without replacement; when a rule resolves to 100% of retained pixels, the metric is calculated once from the full retained-pixel set and repeated across the 50 output rows so the full-pixel condition has zero artificial bootstrap variation. Mean-by-sample-size figures include 95% CI error bars around the 50-iteration mean.
- Source Location: `reports/tables/sample_size_effects/pca_mean_distance/`, `reports/tables/sample_size_effects/spectral_rao_q/`, and `reports/tables/sample_size_effects/alpha_hull_area/`
- Row Count: each metric folder has 392 rows in `*_sample_size_design.csv`, 19,600 rows in `*_sample_size_boot_long.csv`, and 392 rows in `*_sample_size_summary.csv`
- Column Count: 18 columns in design; 20 columns in long bootstrap results; 29 columns in summary
- Primary Keys: `scale`, `quad_id`, `sample_rule`; plus `bootstrap_iter` for long bootstrap results
- Date Fields: none detected
- Numeric Fields: retained pixel counts, sample sizes, sample fractions, bootstrap metric values, bootstrap means/SD/CV/SE, 95 percent confidence intervals, and differences from the fixed-4,000 rule
- Categorical Fields: `scale`, `quad_id`, `sample_rule`, `rule_type`, `rule_axis_label`, `sample_label`, `plot_sample_label`, `metric_method`
- Missing Value Summary: no missing bootstrap mean values in current summary tables
- Last Updated: 2026-07-04

| Variable Group | Type | Description |
|---|---|---|
| `scale`, `quad_id`, `sample_rule`, `rule_type`, `sample_label` | identifier/categorical | Quadrat scale and sample-size rule, with fixed-count labels showing the realized percent of each quadrat in parentheses |
| `n_pixels`, `sample_size`, `sample_fraction` | numeric | Retained pixel pool and number/fraction of pixels sampled per bootstrap iteration |
| `bootstrap_iter`, `metric_value` | numeric | Bootstrap iteration index and resulting PCA mean distance, Rao's Q, or alpha-hull area value, depending on the metric folder |
| `metric_mean`, `metric_sd`, `metric_cv`, `metric_se` | numeric | Mean and variability of the 50 bootstrap metric values for each quadrat x sample-rule combination |
| `ci_low`, `ci_high` | numeric | 95 percent confidence interval around the 50-iteration mean |
| `delta_from_fixed_4000` | numeric | Difference between each sample rule mean and the same quadrat's fixed-4,000-pixel mean |
| `metric_method`, `metric_n_points_min`, `metric_n_points_max` | categorical/numeric | Calculation method and point counts used by the metric; alpha-hull may remove duplicate PC1-PC2 points internally |

### Current all-scale plant diversity outputs

- Dataset Name: Current per-scale species and phylogenetic diversity outputs
- Description: Quadrat-level species matrix, species diversity metrics, and phylogenetic diversity metrics generated by `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`. The species matrix values represent the summed proportion of buffered tree crowns intersecting each quadrat by species. The workflow preserves the current tree filtering logic from the older plant-diversity scripts: `DBH.2024 >= 200`, `crown.position` in `3, 4, 5`, and `cluster_status` in `A, R`.
- Source Location: `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv`, `Quad_Values/Diversity_SHPs/plant_diversity_20m.csv`, `Quad_Values/Diversity_SHPs/plant_diversity_50m.csv`, with matching shapefiles
- Row Count: 2,000 rows for 10 m, 500 rows for 20 m, and 80 rows for 50 m
- Column Count: 61 columns for 10 m CSV, 60 columns for 20 m CSV, and 60 columns for 50 m CSV
- Primary Keys: `scale`, `quad_id`
- Date Fields: none detected
- Numeric Fields: species-code crown-overlap values, `richness`, `shannon`, `simpson`, `evenness`, `faith_pd`, `rao_pd`, `afaith_pd`
- Categorical Fields: `quad_id`, `scale`, and original quadrat ID fields such as `Name` and `sub_id`
- Missing Value Summary: no missing values in required diversity metrics; every current spectral `quad_id` has a matching plant-diversity row
- Last Updated: 2026-06-22

| Variable Group | Type | Description |
|---|---|---|
| `quad_id`, `scale` | identifier | Scale-aware quadrat identifier aligned with current spectral diversity outputs: 10 m `sub_id`, 20 m `Name`, and 50 m `Name` |
| `Name`, `sub_id` | identifier | Original quadrat shapefile identifiers retained where available |
| species-code columns such as `ACRU`, `LITU`, `QURU` | numeric | Summed crown-overlap proportions for each species in each quadrat |
| `richness` | numeric | Count of species with positive crown-overlap values |
| `shannon` | numeric | Shannon species diversity calculated from positive species crown-overlap values |
| `simpson` | numeric | Simpson diversity calculated from positive species crown-overlap values |
| `evenness` | numeric | Shannon diversity divided by log richness when richness is greater than 1 |
| `faith_pd` | numeric | Faith's phylogenetic diversity from the retained species in each quadrat |
| `rao_pd` | numeric | Rao phylogenetic diversity from species crown-overlap weights and phylogenetic distances |
| `afaith_pd` | numeric | Abundance-weighted Faith's phylogenetic diversity using species crown-overlap weights |

### Current all-scale environmental variable outputs

- Dataset Name: Current per-scale elevation and topographic roughness outputs
- Description: Quadrat-level environmental covariates generated by `scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R` from `PRFPD_DTM_leafOff.tiff`. The workflow calculates mean DTM elevation and mean Riley topographic roughness index values after calculating Riley TRI rasters with 5x5, 11x11, and 21x21 moving windows.
- Source Location: `Quad_Values/Enviro_SHPs/enviro_variables_10m.csv`, `enviro_variables_20m.csv`, and `enviro_variables_50m.csv`, with matching shapefiles
- Row Count: 2,000 rows for 10 m, 500 rows for 20 m, and 80 rows for 50 m
- Column Count: 8 columns for 10 m CSV, 7 columns for 20 m CSV, and 7 columns for 50 m CSV
- Primary Keys: `scale`, `quad_id`
- Date Fields: none detected
- Numeric Fields: `elev_mean`, `tri5_mean`, `tri11_mean`, `tri21_mean`, and source quadrat identifier fields such as `Name`
- Categorical Fields: `quad_id`, `scale`, and 10 m `sub_id`
- Missing Value Summary: no missing values in the required environmental metrics
- Last Updated: 2026-06-23

| Variable | Type | Description |
|---|---|---|
| `quad_id` | identifier | Scale-aware quadrat identifier aligned with current biodiversity and spectral outputs: 10 m `sub_id`, 20 m `Name`, and 50 m `Name` |
| `scale` | categorical | Quadrat scale: `10m`, `20m`, or `50m` |
| `elev_mean` | numeric | Mean DTM elevation within the quadrat, in meters |
| `tri5_mean` | numeric | Mean Riley topographic roughness index calculated from a 5x5 DTM moving window |
| `tri11_mean` | numeric | Mean Riley topographic roughness index calculated from an 11x11 DTM moving window |
| `tri21_mean` | numeric | Mean Riley topographic roughness index calculated from a 21x21 DTM moving window |
| `Name`, `sub_id` | identifier | Original quadrat shapefile identifiers retained where available |

### Current combined quadrat analysis tables

- Dataset Name: Current combined per-scale quadrat analysis tables
- Description: Analysis-ready root-level CSVs combining current spectral heterogeneity values, species composition, species diversity, phylogenetic diversity, environmental/topographic values, and quadrat center coordinates. Generated by `scripts/3_Analysis/combine_quadrat_analysis_tables.R`.
- Source Location: `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`
- Row Count: 2,000 rows for 10 m, 500 rows for 20 m, and 80 rows for 50 m
- Column Count: 24 columns for each scale
- Primary Keys: `scale`, `quad_id`
- Date Fields: none detected
- Numeric Fields: `center_x`, `center_y`, species diversity metrics, phylogenetic diversity metrics, spectral heterogeneity metrics, and environmental/topographic metrics
- Categorical Fields: `quad_id`, `scale`
- Missing Value Summary: center coordinate and environmental fields have no missing values; spectral fields retain upstream missingness where current raster summaries are absent or PCA-dependent values were manually excluded upstream. Validation found missing `spec_sa` values: 97 at 10 m, 15 at 20 m, and 0 at 50 m. Missing `spec_pca_mean` values: 256 at 10 m, 64 at 20 m, and 6 at 50 m.
- Last Updated: 2026-06-24

| Variable Group | Type | Description |
|---|---|---|
| `quad_id`, `scale` | identifier | Scale-aware quadrat identifier and quadrat grain |
| `center_x`, `center_y` | numeric/spatial | Quadrat polygon centroid coordinates from the plant-diversity shapefile CRS, NAD83 / UTM zone 16N |
| `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even` | numeric | Species richness, Shannon diversity, Simpson diversity, and evenness |
| `phy_faith`, `phy_rao`, `phy_afaith` | numeric | Faith's PD, Rao PD, and abundance-weighted Faith's PD |
| `spec_sa`, `spec_pca_mean`, `spec_pca_med`, `spec_pca_sd`, `spec_rao_q`, `spec_alpha`, `spec_convex`, `spec_hull3d_v`, `spec_hull3d_a` | numeric | Shortened spectral heterogeneity metrics from current smoothed, masked 5 nm spectral outputs |
| `env_elev`, `env_tri5`, `env_tri11`, `env_tri21` | numeric | Mean DTM elevation and mean Riley TRI values from 5x5, 11x11, and 21x21 windows |

Detailed definitions and calculation notes for every shortened column are in `reports/combined_quadrat_variable_guide.md` and the user-friendly Word version `combined_quadrat_variable_guide.docx`.

### `reports/tables/multiscale_spectral_biodiversity/*.csv`

- Dataset Name: Multiscale spectral-biodiversity analysis outputs
- Description: Derived analysis tables generated by `scripts/3_Analysis/multiscale_spectral_biodiversity_analysis.R` from the current root-level combined quadrat analysis tables. These tables support the PDF findings report and model appendix.
- Source Location: `reports/tables/multiscale_spectral_biodiversity/`
- Row Count: 2,580 rows in `analysis_dataset_with_flags.csv`; 60 rows in `correlation_results.csv`; 105 rows in `model_summary_results.csv`; 309 rows in `model_coefficient_results.csv`; 3 rows in `primary_best_models_by_scale.csv`; 9 rows in `residual_moran_diagnostics.csv`
- Column Count: varies by output table
- Primary Keys: `scale` plus `quad_id` for the flagged analysis dataset; `scale`, `response`, `model_type`, and `predictor` for model summaries; `scale`, `response`, and `predictor` for correlations
- Date Fields: none detected
- Numeric Fields: standardized analysis variables, Pearson/Spearman correlations, model R2 and AIC values, coefficients, confidence intervals, p-values, and Moran's I diagnostics
- Categorical Fields: `quad_id`, `scale`, `response`, `predictor`, `model_type`, edge/primary-analysis flags
- Missing Value Summary: preserves upstream missingness from the combined quadrat tables; no spectral values are imputed
- Last Updated: 2026-06-24

| Variable Group | Type | Description |
|---|---|---|
| `quad_id`, `scale`, `edge_flag`, `primary_analysis` | identifier/categorical | Quadrat identity, analysis scale, documented edge-quadrat flag, and whether the row is included in primary model inference |
| `*_z` fields | numeric | Standardized model variables used for comparable coefficients |
| `pearson_r`, `spearman_r`, `pearson_p` | numeric | Correlation diagnostics by scale, spectral response, and predictor |
| `r_squared`, `adj_r_squared`, `aic`, `delta_aic` | numeric | Candidate-model performance summaries |
| `estimate`, `std.error`, `conf.low`, `conf.high`, `p.value` | numeric | Model coefficient summaries from standardized OLS candidate models |
| `moran_i`, `expected_i` | numeric | 8-nearest-neighbor Moran's I diagnostics for residual spatial autocorrelation |

### `Documents/PDFs/spectral_biodiversity_*.pdf`

- Dataset Name: Multiscale spectral-biodiversity PDF reports
- Description: User-facing PDF reports summarizing the current multiscale spectral-biodiversity analysis and model appendix.
- Source Location: `Documents/PDFs/spectral_biodiversity_multiscale_findings.pdf` and `Documents/PDFs/spectral_biodiversity_model_appendix.pdf`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: not applicable
- Date Fields: file timestamps only
- Numeric Fields: embedded model, correlation, and diagnostic summaries
- Categorical Fields: report sections, captions, and interpretation notes
- Missing Value Summary: not applicable
- Last Updated: 2026-06-24

### `Quad_Values/Other_variables/Slope_Aspect.csv`

- Dataset Name: 20 m slope/aspect environmental covariates
- Description: Quadrat-level topographic covariates used as environmental drivers.
- Source Location: `Quad_Values/Other_variables/Slope_Aspect.csv`
- Row Count: 500
- Column Count: 7
- Primary Keys: `name`
- Date Fields: none detected
- Numeric Fields: unnamed first column, `name`, `slope_percent`, `aspect_sin`, `aspect_cos`, `tei_north`, `tei_east`
- Categorical Fields: none detected
- Missing Value Summary: none detected
- Last Updated: 2026-05-05

| Variable | Type | Description |
|---|---|---|
| unnamed first column | numeric | Likely row/index column |
| `name` | numeric | Quadrat identifier |
| `slope_percent` | numeric | Slope percentage |
| `aspect_sin` | numeric | Sine-transformed aspect |
| `aspect_cos` | numeric | Cosine-transformed aspect |
| `tei_north` | numeric | Topographic exposure index, north component |
| `tei_east` | numeric | Topographic exposure index, east component |

## Spatial and Raster Dataset Groups

### `HSI_NA_trimmed/`

- Dataset Name: Trimmed hyperspectral imagery
- Description: Current raw/trimmed hyperspectral raster cubes used as primary spectral inputs.
- Source Location: `HSI_NA_trimmed/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: raster tile names
- Date Fields: file timestamps only
- Numeric Fields: raster spectral bands
- Categorical Fields: none
- Missing Value Summary: not inspected; requires geospatial raster tooling
- Last Updated: 2026-03-02

### `Quad_Scale_SHPs/`

- Dataset Name: Quadrat scale shapefiles
- Description: Spatial partitioning layers for 10 m, 20 m, and 50 m quadrats.
- Source Location: `Quad_Scale_SHPs/`
- Row Count: not inspected from shapefile DBF in this pass
- Column Count: not inspected from shapefile DBF in this pass
- Primary Keys: quadrat identifiers in shapefile attributes
- Date Fields: none expected
- Numeric Fields: quadrat IDs and geometry attributes
- Categorical Fields: geometry metadata
- Missing Value Summary: not inspected
- Last Updated: 2026-01-27

### `Quad_Spectra/`

- Dataset Name: Quadrat-level spectral extracts
- Description: Current quadrat-level spectral rasters for 10 m, 20 m, and 50 m scales. The user has independently tested the partitioned quadrat spectra and confirmed that `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m` are the spectral inputs to use moving forward. Current Savitzky-Golay smoothed outputs derived from those inputs are stored in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`. Current 5 nm resampled outputs derived from those smoothed spectra are stored in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`. Additional folders contain older smoothed products, RGB products, vegetation-index products, and testing artifacts.
- Source Location: `Quad_Spectra/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: quadrat raster filenames
- Date Fields: file timestamps only
- Numeric Fields: raster bands and derived index values
- Categorical Fields: none
- Missing Value Summary: not inspected; requires raster reads
- Last Updated: 2026-06-18

### `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, `Quad_Spectra/50m_smooth`

- Dataset Name: Current smoothed quadrat spectra
- Description: Savitzky-Golay smoothed ENVI raster outputs generated from the confirmed partitioned spectra in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Source Location: `Quad_Spectra/*_smooth`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: quadrat raster filenames
- Date Fields: file timestamps only
- Numeric Fields: 328 spectral bands per validated output raster
- Categorical Fields: none
- Missing Value Summary: not summarized globally; sample outputs open successfully with `terra::rast()`
- Last Updated: 2026-06-17

| Scale | Output Folder | Raster Count | Approx. Folder Size | Validation Sample |
|---|---|---:|---:|---|
| 10 m | `Quad_Spectra/10m_smooth` | 1,909 | 46.564 GB | `0_a`: 328 bands, 136 rows, 136 columns |
| 20 m | `Quad_Spectra/20m_smooth` | 485 | 47.133 GB | `0`: 328 bands, 272 rows, 272 columns |
| 50 m | `Quad_Spectra/50m_smooth` | 80 | 48.544 GB | `sub50_1`: 328 bands, 680 rows, 680 columns |

### `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, `Quad_Spectra/50m_smooth_5nm`

- Dataset Name: Current smoothed 5 nm quadrat spectra
- Description: ENVI raster outputs generated by resampling current smoothed spectra to 5 nm spectral spacing.
- Source Location: `Quad_Spectra/*_smooth_5nm`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: quadrat raster filenames
- Date Fields: file timestamps only
- Numeric Fields: 121 spectral bands per validated output raster, from 398 nm through 998 nm
- Categorical Fields: none
- Missing Value Summary: not summarized globally; sample outputs open successfully with `terra::rast()`
- Last Updated: 2026-06-18

| Scale | Output Folder | Raster Count | Approx. Folder Size | Validation Sample |
|---|---|---:|---:|---|
| 10 m | `Quad_Spectra/10m_smooth_5nm` | 1,909 | 17.177 GB | `0_a`: 121 bands, 136 rows, 136 columns |
| 20 m | `Quad_Spectra/20m_smooth_5nm` | 485 | 17.387 GB | `0`: 121 bands, 272 rows, 272 columns |
| 50 m | `Quad_Spectra/50m_smooth_5nm` | 80 | 17.908 GB | `sub50_1`: 121 bands, 680 rows, 680 columns |

### `Quad_Values/Diversity_Rasters/`

- Dataset Name: Diversity raster outputs
- Description: Rasterized biodiversity metrics, including Shannon, Simpson, richness, Rao PD, Faith PD, abundance-weighted Faith PD, and evenness products.
- Source Location: `Quad_Values/Diversity_Rasters/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: raster filenames
- Date Fields: file timestamps only
- Numeric Fields: raster cell values
- Categorical Fields: none
- Missing Value Summary: not inspected
- Last Updated: based on file metadata

### `Quad_Values/Diversity_SHPs/`

- Dataset Name: Diversity shapefile outputs
- Description: Vector biodiversity metrics for analysis quadrats.
- Source Location: `Quad_Values/Diversity_SHPs/`
- Row Count: not inspected from shapefile DBF in this pass
- Column Count: not inspected from shapefile DBF in this pass
- Primary Keys: quadrat identifiers
- Date Fields: none expected
- Numeric Fields: biodiversity metric fields
- Categorical Fields: geometry metadata
- Missing Value Summary: not inspected
- Last Updated: based on file metadata

### `Quad_Values/Spectral_diversitySHPs/`

- Dataset Name: Spectral diversity shapefile outputs
- Description: Current vector spectral products that produce spectral heterogeneity values, including spectral angle entropy, band entropy, quadrat-inclusive entropy variants, and global PCA metrics.
- Source Location: `Quad_Values/Spectral_diversitySHPs/`
- Row Count: not inspected from shapefile DBF in this pass
- Column Count: not inspected from shapefile DBF in this pass
- Primary Keys: quadrat identifiers
- Date Fields: none expected
- Numeric Fields: spectral diversity metric fields
- Categorical Fields: geometry metadata
- Missing Value Summary: not inspected
- Last Updated: based on file metadata

### `VegIndex_NA_trimmed/`

- Dataset Name: Trimmed vegetation-index rasters
- Description: Vegetation-index products generated from trimmed hyperspectral imagery.
- Source Location: `VegIndex_NA_trimmed/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: raster tile names
- Date Fields: file timestamps only
- Numeric Fields: vegetation index raster values
- Categorical Fields: none
- Missing Value Summary: not inspected
- Last Updated: 2026-03-03

## Notes

- PowerShell assigned `H1` to unnamed first columns in several CSVs. These are documented as unnamed row/index columns and should not be interpreted as analytical variables unless script review proves otherwise.
- Spectral products that produce spectral heterogeneity values are stored in `Quad_Values/Spectral_diversitySHPs/`.
- Confirmed spectra cropped to quadrat extents are in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current smoothed spectra generated from confirmed partitioned inputs are in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`.
- Current smoothed 5 nm spectra are in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- The next spectral heterogeneity workflow should consume the current `_smooth_5nm` spectra. Its shadow mask retains pixels greater than `0.0305476` at `563 nm` under the current direction setting.
- `Quad_Spectra/10m_test`, `Quad_Spectra/20m_test`, and `Quad_Spectra/50m_test` are testing/validation artifacts unless the user explicitly promotes them.
- Raster, ENVI, and shapefile attribute-level summaries should be expanded later using R packages such as `terra`, `sf`, and `dplyr` once R is available on the execution path.
- Authorized `old/`, `Outdated/`, and `Currently not relevant/` directories were removed on 2026-06-15.
