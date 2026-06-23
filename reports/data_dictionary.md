# Data Dictionary

Last updated: 2026-06-22

This dictionary inventories current datasets visible in the repository. Large raster and hyperspectral products were documented from file metadata and project context; row/column counts are listed where lightweight tabular inspection was possible. Unnamed first columns in CSV outputs are treated as likely index columns unless later script review shows otherwise.

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

### `Indices_SHPs/20m_spectral_sp.csv`

- Dataset Name: 20 m spectral, biodiversity, species, and environmental table
- Description: Main integrated 20 m analysis table combining spectral entropy/PCA metrics, phylogenetic and species-diversity metrics, species abundance columns, geometry, and environmental variables.
- Source Location: `Indices_SHPs/20m_spectral_sp.csv`
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

### `Indices_SHPs/20m_SA_smooth_masked.csv`

- Dataset Name: 20 m masked smoothed spectral angle entropy
- Description: Quadrat-level spectral entropy table using smoothed and masked spectra.
- Source Location: `Indices_SHPs/20m_SA_smooth_masked.csv`
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

### `Indices_SHPs/20m_SA_smooth_masked_7_11.csv`

- Dataset Name: 20 m masked smoothed spectral angle entropy, 7/11 variant
- Description: Quadrat-level spectral entropy table using a smoothed/masked parameter variant.
- Source Location: `Indices_SHPs/20m_SA_smooth_masked_7_11.csv`
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

### `Indices_SHPs/20m_SA_entrop_boot_results.csv`

- Dataset Name: 20 m spectral angle entropy bootstrap results, 10 iterations
- Description: Quadrat-level bootstrap estimates for spectral angle entropy.
- Source Location: `Indices_SHPs/20m_SA_entrop_boot_results.csv`
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
| `boot_1` to `boot_10` | numeric | Bootstrap replicate entropy estimates |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary statistics |

### `Indices_SHPs/20m_SA_entrop_boot100_results.csv`

- Dataset Name: 20 m spectral angle entropy bootstrap results, 100 iterations
- Description: Current expanded bootstrap estimates for spectral angle entropy.
- Source Location: `Indices_SHPs/20m_SA_entrop_boot100_results.csv`
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
| `boot_1` to `boot_100` | numeric | Bootstrap replicate entropy estimates |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary statistics |

### Current spectral angle entropy outputs from `_smooth_5nm`

- Dataset Name: Current per-scale spectral angle entropy from smoothed 5 nm spectra
- Description: Per-quadrat spectral heterogeneity tables generated by `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R` from current smoothed 5 nm spectra. The workflow attempts exact all-pixel spectral angle entropy only when pairwise angle counts are below the configured limit and otherwise uses the mean of bootstrap subsamples as the primary spectral heterogeneity value. Large bootstrap replicates use sampled pixel pairs to avoid materializing all pairwise angles.
- Source Location: `Indices_SHPs/*_SA_entropy_smooth_masked_5nm_*.csv`
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
| `n_pairs` | numeric | Number of possible pairwise spectral angles among retained pixels |
| `method` | categorical | `exact_all_pixels`, `bootstrap_mean`, or `insufficient_pixels` |
| `spectral_entropy` | numeric | Primary spectral heterogeneity value; exact entropy when feasible, otherwise bootstrap mean |
| `exact_entropy` | numeric | All-pixel spectral angle entropy when pair counts are below the configured exact threshold |
| `boot_mean`, `boot_sd`, `boot_cv`, `boot_median`, `boot_min`, `boot_max` | numeric | Bootstrap summary fields used when all-pixel entropy is not computationally reasonable |
| `boot_1` ... `boot_n` | numeric | Wide-format bootstrap replicate values, where `n` is controlled by `SA_N_BOOT` |

### Current combined spectral heterogeneity outputs from `_smooth_5nm`

- Dataset Name: Current per-scale combined spectral heterogeneity outputs
- Description: Exclusion-aware spectral heterogeneity tables generated by `scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R`. These combine existing SA entropy with PCA-dependent metrics. The current PCA basis excludes manually documented atmospheric/cloud-affected quadrats before PCA sampling and leaves PCA-dependent values missing for those quadrats. Earlier PCA-dependent outputs from 2026-06-22 created before this exclusion policy should be disregarded.
- Source Location: `Indices_SHPs/*_spectral_heterogeneity_smooth_masked_5nm_summary.csv` and `Indices_SHPs/Spectral_diversitySHPs/spectral_heterogeneity_*_smooth_masked_5nm.*`
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
| `sa_*` fields | numeric/categorical | Prefixed fields imported from the current SA entropy summary CSVs |

### `reports/tables/bootstrap_variation/*.csv`

- Dataset Name: Bootstrap variation quality-control tables
- Description: Derived validation tables comparing within-quadrat bootstrap variation with between-quadrat spectral heterogeneity variation for 10 m, 20 m, and 50 m spectral entropy outputs.
- Source Location: `reports/tables/bootstrap_variation/`
- Row Count: 3 rows in `bootstrap_variation_scale_summary.csv`; 2,422 rows in `bootstrap_variation_quadrat_diagnostics.csv`; 45 rows in `bootstrap_variation_top_unstable_quadrats.csv`
- Column Count: 36 columns in scale summary; 18 columns in quadrat diagnostics; 18 columns in top unstable quadrats
- Primary Keys: `scale` for scale summary; `scale` + `quad_id` for quadrat-level diagnostics
- Date Fields: none detected
- Numeric Fields: bootstrap SD/CV, between-scale variance summaries, standard-error estimates, pixel counts, pair counts, outlier replicate counts
- Categorical Fields: `scale`, `quad_id`, `method`, `stability_class`
- Missing Value Summary: exact all-pixel and insufficient-pixel quadrats are excluded from quadrat-level bootstrap diagnostics; no missing values expected in the scale summary
- Last Updated: 2026-06-18

| Variable Group | Type | Description |
|---|---|---|
| `scale`, `quad_id`, `method` | identifier/categorical | Analysis scale, quadrat identifier, and spectral entropy estimation method |
| `between_sd`, `between_var`, `between_iqr` | numeric | Between-quadrat spectral entropy variation at each scale |
| `mean_entropy`, `median_entropy`, `entropy_skewness` | numeric | Scale-level distribution summaries for final spectral heterogeneity values |
| `within_sd_*`, `within_cv_*` | numeric | Within-quadrat bootstrap variation summaries |
| `boot_se_*` | numeric | Standard error summaries for the reported bootstrap mean |
| `within_to_between_var_ratio`, `icc_like_between_fraction` | numeric | Diagnostics comparing within-bootstrap variation to between-quadrat variation |
| `quads_cv_gt_005`, `quads_cv_gt_010` | numeric | Counts of quadrats exceeding 5% and 10% bootstrap CV thresholds |
| `outlier_reps` | numeric | Number of bootstrap replicates outside the quadrat-level 1.5 IQR fences |

### Current all-scale plant diversity outputs

- Dataset Name: Current per-scale species and phylogenetic diversity outputs
- Description: Quadrat-level species matrix, species diversity metrics, and phylogenetic diversity metrics generated by `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`. The species matrix values represent the summed proportion of buffered tree crowns intersecting each quadrat by species. The workflow preserves the current tree filtering logic from the older plant-diversity scripts: `DBH.2024 >= 200`, `crown.position` in `3, 4, 5`, and `cluster_status` in `A, R`.
- Source Location: `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.csv`, `Indices_SHPs/Diversity_SHPs/plant_diversity_20m.csv`, `Indices_SHPs/Diversity_SHPs/plant_diversity_50m.csv`, with matching shapefiles
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

### `Indices_SHPs/Other_variables/Slope_Aspect.csv`

- Dataset Name: 20 m slope/aspect environmental covariates
- Description: Quadrat-level topographic covariates used as environmental drivers.
- Source Location: `Indices_SHPs/Other_variables/Slope_Aspect.csv`
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

### `Indices_SHPs/Diversity_Rasters/`

- Dataset Name: Diversity raster outputs
- Description: Rasterized biodiversity metrics, including Shannon, Simpson, richness, Rao PD, Faith PD, abundance-weighted Faith PD, and evenness products.
- Source Location: `Indices_SHPs/Diversity_Rasters/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: raster filenames
- Date Fields: file timestamps only
- Numeric Fields: raster cell values
- Categorical Fields: none
- Missing Value Summary: not inspected
- Last Updated: based on file metadata

### `Indices_SHPs/Diversity_SHPs/`

- Dataset Name: Diversity shapefile outputs
- Description: Vector biodiversity metrics for analysis quadrats.
- Source Location: `Indices_SHPs/Diversity_SHPs/`
- Row Count: not inspected from shapefile DBF in this pass
- Column Count: not inspected from shapefile DBF in this pass
- Primary Keys: quadrat identifiers
- Date Fields: none expected
- Numeric Fields: biodiversity metric fields
- Categorical Fields: geometry metadata
- Missing Value Summary: not inspected
- Last Updated: based on file metadata

### `Indices_SHPs/Spectral_diversitySHPs/`

- Dataset Name: Spectral diversity shapefile outputs
- Description: Current vector spectral products that produce spectral heterogeneity values, including spectral angle entropy, band entropy, quadrat-inclusive entropy variants, and global PCA metrics.
- Source Location: `Indices_SHPs/Spectral_diversitySHPs/`
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
- Spectral products that produce spectral heterogeneity values are in `Indices_SHPs/Spectral_diversitySHPs/`.
- Confirmed spectra cropped to quadrat extents are in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`.
- Current smoothed spectra generated from confirmed partitioned inputs are in `Quad_Spectra/10m_smooth`, `Quad_Spectra/20m_smooth`, and `Quad_Spectra/50m_smooth`.
- Current smoothed 5 nm spectra are in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`.
- The next spectral heterogeneity workflow should consume the current `_smooth_5nm` spectra. Its shadow mask retains pixels greater than `0.0305476` at `563 nm` under the current direction setting.
- `Quad_Spectra/10m_test`, `Quad_Spectra/20m_test`, and `Quad_Spectra/50m_test` are testing/validation artifacts unless the user explicitly promotes them.
- Raster, ENVI, and shapefile attribute-level summaries should be expanded later using R packages such as `terra`, `sf`, and `dplyr` once R is available on the execution path.
- Authorized `old/`, `Outdated/`, and `Currently not relevant/` directories were removed on 2026-06-15.
