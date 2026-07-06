# Research Objective and Analysis Specification for Codex Agent

## Project Title

**Evaluating Relationships Between Hyperspectral Spectral Heterogeneity and Phylogenetic Diversity in the Paint Rock Forest Dynamics Plot**

## Purpose

The purpose of this project is to quantify and evaluate relationships between remotely sensed spectral heterogeneity and field-derived biodiversity metrics within the Paint Rock Forest Dynamics Plot. Specifically, this research investigates whether variation in hyperspectral reflectance can serve as a proxy for phylogenetic diversity across multiple spatial scales.

The Codex agent will assist in:

- Reviewing and understanding existing R scripts.
- Identifying dependencies, inputs, and outputs within the workflow.
- Reproducing analyses in a transparent and documented manner.
- Generating summary statistics, figures, tables, and model outputs.
- Evaluating relationships between spectral heterogeneity metrics and biodiversity metrics.
- Documenting assumptions, limitations, and potential sources of error.
- Producing publication-ready outputs where possible.

## Scientific Objectives

### Objective 1: Quantify Spectral Heterogeneity

Calculate spectral heterogeneity metrics from drone-acquired hyperspectral imagery using quadrat-based sampling approaches.

**Primary metric:**

- Spectral entropy derived from pairwise spectral angle distributions. Current SA entropy outputs are usually the mean of 70 bootstrap iterations using up to 5,000 retained pixels per iteration; only the small subset of quadrats below the configured pair-count threshold uses exact all-retained-pixel pairwise angles.

**Potential secondary metrics:**

- Alpha hull of pixels in PCA space
- Euclidean distance in PCA space
- Rao's Quadratic Entropy
- 

**Questions:**

- How does spectral heterogeneity vary across the study area?
- How sensitive are estimates to spatial grain size and pixel selection procedures?

### Objective 2: Quantify Biodiversity Metrics

Generate spatially explicit biodiversity layers derived from forest census data.

**Metrics include:**

- Shannon diversity
- Faith's phylogenetic diversity
- Abundance-weighted Faith's phylogenetic diversity

**Questions:**

- How do biodiversity patterns vary spatially?
- How do diversity estimates change across analysis scales?

### Objective 3: Evaluate Spectral-Biodiversity Relationships

Determine whether spectral heterogeneity is associated with biodiversity metrics.

**Primary hypotheses:**

- H1: Spectral heterogeneity increases with phylogenetic diversity.
- H2: Relationships vary across spatial scales.
- H3: Environmental variables explain additional variation beyond biodiversity metrics alone.

**Expected analyses:**

- Correlation analyses
- Linear regression
- Generalized linear models where appropriate
- Spatial diagnostics
- Model comparison and variance partitioning

**Primary response variables:**

- Spectral entropy and other spectral heterogeneity metrics

**Primary predictor variables:**

- Faith's phylogenetic diversity
- Abundance-weighted Faith's phylogenetic diversity
- Shannon diversity

### Objective 4: Assess Environmental Drivers

Investigate whether environmental variables contribute to unexplained variation in spectral heterogeneity.

**Potential variables:**

- Elevation
- Topographic position
- Additional raster products available for the study area

**Questions:**

- How much variation is explained by environmental factors?
- Does phylogenetic diversity remain significant after environmental controls are included?

### Objective 5: Evaluate Scale Dependence

Assess whether relationships change with spatial grain.

**Current analysis scales:**

- 10 m quadrats
- 20 m quadrats
- 50 m quadrats

**Questions:**

- Which scale produces the strongest biodiversity-spectral relationship?
- Are observed relationships consistent across scales?

## Data Sources

### Data Organization and Processing Structure

#### Raw Spectral Data Location

The raw hyperspectral imagery used in this analysis is stored in:

- **HSI_NA_trimmed/**

This directory contains raw hyperspectral raster files (ENVI format), including paired `.hdr` metadata files and associated auxiliary files. Example entries include:

- raw_0_rd_rf_or
- raw_11159_rd_rf_or
- raw_11892_rd_rf_or

Each file represents a spatially referenced hyperspectral cube used as the primary input for downstream spectral processing.

---

#### Spatial Partitioning Inputs (Quadrat Shapefiles)

Spatial partitioning of hyperspectral data into analysis units is defined in:

- **Quad_Scale_SHPs/**

This directory contains shapefiles used to segment imagery into fixed spatial grains:

- 10 m quadrats (PR_10m.shp)
- 20 m quadrats (PR_20m.shp)
- 50 m quadrats (PR_50m.shp)

Each scale includes standard GIS components (.shp, .dbf, .prj, .shx) and optional KML exports.

These shapefiles define the spatial sampling framework for all spectral and biodiversity aggregation procedures.

---

#### Quadrat-Based Spectral Processing Outputs

Processed and partitioned spectral data are stored in:

- **Quad_Spectra/**

This directory contains quadrat-level spectral extractions and derived products organized by spatial scale (10 m, 20 m, 50 m). The user has independently tested the partitioned quadrat spectra and confirmed that the primary spectra to use moving forward are:

- **Quad_Spectra/10m/**
- **Quad_Spectra/20m/**
- **Quad_Spectra/50m/**

These confirmed folders should be treated as the current spectral inputs for downstream spectral heterogeneity calculations. Testing folders such as **10m_test/**, **20m_test/**, and **50m_test/** are validation artifacts unless the user explicitly promotes them.

Quad_Spectra/ also includes multiple processing stages:

- Raw quadrat spectra (10m/, 20m/, 50m/)
- Current smoothed quadrat spectra (10m_smooth/, 20m_smooth/, 50m_smooth/)
- Current smoothed and 5 nm resampled quadrat spectra (10m_smooth_5nm/, 20m_smooth_5nm/, 50m_smooth_5nm/)
- Spectral smoothing outputs (*_smoothed_5nm, 50m_RGB_smooth)
- Resampled spectra (*_resampled_5nm)
- RGB composites (20m_RGB, 50m_RGB)
- Vegetation indices (*_VegIndex)

These outputs represent the primary analytical dataset used for spectral heterogeneity calculations.

---

#### Quadrat Quality Flags and Analysis Exclusions

Some quadrats require special handling because they represent plot edges or atmospheric distortion artifacts. These flags should be treated as current analysis assumptions unless a later validated script or user instruction supersedes them.

**Edge quadrats**

Edge quadrats are identified in:

- **scripts/3_Analysis/Analysis_PDF.R**

The current 20 m edge quadrat IDs are:

- 0, 1, 100, 2, 200, 3, 300, 4, 400, 5, 500, 6, 600, 7, 700, 8, 800, 9, 900, 10, 1000, 11, 1100, 12, 1200, 13, 14, 1300, 15, 1400, 16, 1500, 17, 1600, 18, 1700, 19, 1800, 20, 1900, 21, 22, 1901, 23, 24, 1903, 124, 1905, 224, 1906, 324, 1907, 424, 1908, 524, 1909, 624, 1910, 1911, 724, 1912, 824, 1913, 924, 1914, 1024, 1915, 1124, 1916, 1224, 1324, 1917, 1424, 1919, 1524, 1920, 1624, 1921, 1724, 1922, 1824, 1923, 1924, 1904, 1902, 1918

These IDs are currently filtered from the 20 m spectral-biodiversity analysis in `Analysis_PDF.R`.

**Atmospheric/cloud-affected quadrats**

Quadrats affected by atmospheric distortion or clouds are identified by the exclusion vector in:

- **scripts/2_Indices Creation/Spectral_diversity/HSI_global_PCA.R**

The current 20 m atmospheric/cloud-affected quadrat IDs are:

- 1424, 1423, 1422, 1420, 1421, 1419, 1418, 1414, 1521, 1522, 1523, 1524, 1520, 1519, 1624, 1622, 1623, 1621, 1620, 1724, 1723, 1722, 1721, 1824, 1823, 1822, 1923, 1924, 1922, 1921, 1322, 1321, 1319, 1320, 1318, 1317, 1316, 1315, 1314, 1313, 1221, 1220, 1219, 1216, 1215, 1213, 1214, 1212, 1211, 1120, 1119, 1115, 1114, 1113, 1112, 1111, 1110, 1014, 1013, 1010, 1009, 909, 908, 24

These IDs are excluded from the 20 m global PCA spectral heterogeneity workflow in `HSI_global_PCA.R`.

**10 m scale rule**

The 10 m quadrats inherit these flags from their parent 20 m quadrats. For each listed 20 m quadrat ID, the corresponding 10 m flagged quadrats are the four suffixed subquadrats:

- `<20m_id>_a`
- `<20m_id>_b`
- `<20m_id>_c`
- `<20m_id>_d`

For example, if 20 m quadrat `1424` is cloud-affected, then `1424_a`, `1424_b`, `1424_c`, and `1424_d` should be treated as cloud-affected at the 10 m scale. These 20 m-derived flags should not be automatically applied to 50 m quadrats unless a separate 50 m-specific rule is documented.

**50 m scale rule**

The current 50 m atmospheric/cloud-affected quadrat IDs are:

- `sub50_80`, `sub50_79`, `sub50_71`, `sub50_70`, `sub50_62`, `sub50_53`

These 50 m quadrats should be treated as manually excluded for PCA-dependent spectral heterogeneity products. They should not contribute pixels to the global PCA basis, and their PCA-dependent heterogeneity values should be left missing in the PCA metric outputs. Spectral angle entropy outputs may still retain their existing values unless a separate SA entropy exclusion workflow is requested.

---

#### Quadrat Extraction Workflow (Scripts)

A scripted workflow is used to partition raw hyperspectral imagery into quadrats using the shapefiles above. These scripts are located in:

- **scripts/**

Key processing stages include:

- **1_Data Processing/** — partitions raw hyperspectral rasters into quadrat-based subsets using Quad_Scale_SHPs
- **2_Indices Creation/** — generates vegetation indices and derived spectral products
- **3_Analysis/** — performs statistical analysis linking spectral metrics to biodiversity
- **auxilary/** — helper functions and utilities
- **Outdated/** — deprecated or legacy workflows
- **visuals/** — figure generation scripts

The quadrat partitioning workflow produces the organized outputs stored in Quad_Spectra/, which serve as the input for downstream spectral heterogeneity and biodiversity analyses.

###Script where quad partioning occurs: /scripts/1_Data Processing/HSI_quad_crop_refined.R

---

### Forest Census Data

Contains:

- Species identity
- Tree coordinates
- DBH
- Canopy Position Index (CPI)
- Tag identifiers

**Filtered analysis population:**

- CPI = 4 or 5 OR
- DBH > 200 mm

**Expected dataset size:**

- Approximately 4,519 upper canopy individuals
- Approximately 51 species

### Hyperspectral Imagery

**Sensor:**

- Headwall Micro A-Series VNIR

**Characteristics:**

- 396–1002 nm spectral range
- 328 bands
- Approximately 2 nm spectral resolution
- 10 cm spatial resolution

**Derived products:**

- Reflectance mosaics
- Shadow masks
- Spectral heterogeneity rasters

## Codex Agent Responsibilities

### Script Review

For each R script:

- Identify purpose.
- List required packages.
- Document inputs.
- Document outputs.
- Identify assumptions.
- Identify computational bottlenecks.
- Detect unused code or variables.
- Include identification of whether inputs originate from raw hyperspectral data (HSI_NA_trimmed/) or pre-partitioned quadrat outputs (Quad_Spectra/).

### Reproducibility

The agent should:

- Prefer modular functions.
- Review existing scripts in **scripts/** before writing new R code, especially loop structures used for iterating across quadrats, scales, and spectral files.
- Avoid hard-coded file paths.
- Clearly document workflow dependencies.
- Record package versions when possible.
- Flag analyses that cannot be reproduced due to missing data.

### Output Generation

Produce:

- Data tables (.CSV)
- Summary tables
- Regression summaries
- Diagnostic plots
- Spatial maps
- Correlation matrices
- Variable importance summaries

All outputs should be saved using consistent naming conventions.

## Expected Deliverables

### Figures

Potential figures include:

- Study area map
- Spectral heterogeneity maps
- Biodiversity maps
- Scale comparison figures
- Spectral diversity versus phylogenetic diversity scatterplots
- Environmental covariate maps
- Model diagnostic plots

### Tables

Potential tables include:

- Species composition summary
- Diversity metric summary statistics
- Correlation matrix
- Model coefficient tables
- Variance explained by predictor groups

### Statistical Outputs

Report:

- Coefficients
- Standard errors
- Confidence intervals
- R² values
- Adjusted R² values
- P-values
- Model diagnostics

## Success Criteria

The analysis will be considered successful if:

- All scripts execute reproducibly.
- Spectral heterogeneity metrics are generated correctly.
- Biodiversity metrics are reproduced from census data.
- Relationships between spectral and biodiversity metrics are quantified.
- Scale-dependent effects are evaluated.
- Environmental influences are assessed.
- Results are documented sufficiently for manuscript preparation and future replication.

## Notes for Agent

This project is exploratory but hypothesis-driven. Preference should be given to transparent, reproducible workflows over computational complexity. Any modifications to existing analytical methods should be clearly documented and justified. The agent should maintain a detailed record of assumptions, data transformations, filtering procedures, and statistical decisions throughout the workflow.
