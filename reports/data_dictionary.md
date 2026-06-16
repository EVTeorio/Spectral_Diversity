# Data Dictionary

Last updated: 2026-06-15

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
- Description: Current quadrat-level spectral rasters for 10 m, 20 m, and 50 m scales, including raw subsets, smoothed products, resampled products, RGB products, and vegetation-index products.
- Source Location: `Quad_Spectra/`
- Row Count: not applicable
- Column Count: not applicable
- Primary Keys: quadrat raster filenames
- Date Fields: file timestamps only
- Numeric Fields: raster bands and derived index values
- Categorical Fields: none
- Missing Value Summary: not inspected; requires raster reads
- Last Updated: 2026-04-01

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
- Spectra cropped to quadrat extents are in `Quad_Spectra/`.
- Raster, ENVI, and shapefile attribute-level summaries should be expanded later using R packages such as `terra`, `sf`, and `dplyr` once R is available on the execution path.
- Authorized `old/`, `Outdated/`, and `Currently not relevant/` directories were removed on 2026-06-15.
