# Directory Map

Last updated: 2026-06-15

Working directory:
`C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity`

## Top-Level Inventory

| Directory | Purpose | File Count | Approx. Size |
|---|---|---:|---:|
| `Documents/` | Manuscript drafts, reports, presentations, and workshop materials | 18 | 0.015 GB |
| `HSI_NA_trimmed/` | Current trimmed raw hyperspectral imagery in ENVI-style raster format | 75 | 57.022 GB |
| `Indices_SHPs/` | Current derived biodiversity, spectral diversity, and environmental index outputs | 71 | 0.024 GB |
| `logs/` | Project logs required by governance standards | 0 | 0 GB |
| `Quad_Scale_SHPs/` | Quadrat boundary shapefiles and KMLs for 10 m, 20 m, and 50 m analysis scales | 14 | 0.002 GB |
| `Quad_Spectra/` | Quadrat-level spectral raster extracts and derived products | 30,035 | 252.771 GB |
| `reports/` | Governance reports, project inventories, task reports, plans, and validation notes | 7 | <0.001 GB |
| `scripts/` | R workflow scripts for preprocessing, index creation, analysis, utilities, and visuals | 47 | <0.001 GB |
| `Shadow_vs_sunlit_SHP/` | Shadow/sunlit training or classification shapefiles | 12 | <0.001 GB |
| `Test/` | Ad hoc test/projection script area | 1 | <0.001 GB |
| `tests/` | Formal unit tests for future `testthat` coverage using representative data subsamples | 0 | 0 GB |
| `VegIndex_NA_trimmed/` | Vegetation-index raster outputs aligned with trimmed imagery | 75 | 16.485 GB |
| `Visuals/` | Current exported figures and visual raster products | 5 | 0.001 GB |

## Repository Structure

```text
Spectral_Diversity/
    CODEX_AGENT_GUIDELINES.md
    RESEARCH_OBJECTIVES.md
    markdown_script.R
    global_pca.rds
    PR_tree_DL.csv
    51sp_taxanomy.csv
    species_summary.csv
    spec_shannaon_DIV_with Geometry.csv
    51sp_codeTree.pdf
    51sp_codeTree.png
    PlantTrait_Phylogene_Chart.drawio
    hsi_illumination_class.tif.aux.xml

    Documents/
        manuscript drafts, PDF model reports, ARD presentation, workshop document

    HSI_NA_trimmed/
        trimmed hyperspectral raster cubes and paired metadata

    Indices_SHPs/
        20m_spectral_sp.csv
        20m_SA_smooth_masked.csv
        20m_SA_smooth_masked_7_11.csv
        20m_SA_entrop_boot_results.csv
        20m_SA_entrop_boot100_results.csv
        Diversity_Rasters/
        Diversity_SHPs/
        Other_variables/
            Slope_Aspect.csv
            Slope_Aspect_rasters/
            lci_update/
        Spectral_diversitySHPs/

    Quad_Scale_SHPs/
        PR_10m.*
        PR_20m.*
        PR_50m.*
        KMLs/

    Quad_Spectra/
        10m/
        10m_resampled_5nm/
        10m_smoothed_5nm/
        10m_VegIndex/
        20m/
        20m_resampled_5nm/
        20m_RGB/
        20m_smoothed_5nm/
        20m_VegIndex/
        50m/
        50m_resampled_5nm/
        50m_RGB/
        50m_RGB_smooth/
        50m_smoothed_5nm/
        50m_VegIndex/

    scripts/
        1_Data Processing/
        2_Indices Creation/
            Plant_diversity/
            Spectral_diversity/
        3_Analysis/
        auxilary/
        visuals/

    Shadow_vs_sunlit_SHP/
        Shadow_vs_sunlit.*
        Shadow_vs_sunlit_Rev.*

    Test/
        UTM projected.R

    tests/

    VegIndex_NA_trimmed/
        vegetation-index raster outputs

    Visuals/
        exported figures and map products

    reports/
        architecture/
        execution_plans/
        history/
        tasks/
        validation/

    logs/
```

## Current Script Inventory

### `scripts/1_Data Processing/`

- `aspect_slope.R`
- `HSI_NA_trim.R`
- `HSI_quad_crop_refined.R`
- `resampling.R`
- `smoothing.R`
- `Sunlit_vs_Shadow_AUC.R`

### `scripts/2_Indices Creation/`

- `adding_other_variables.R`
- `adding_to_CSV.R`
- `elevation.R`
- `Indices_CSV.R`
- `topo_Index.R`
- `Plant_diversity/Missing_quads.R`
- `Plant_diversity/phylogenetic_diversity.R`
- `Plant_diversity/sp_weighted_matrix.R`
- `Plant_diversity/species_diversity.R`
- `Spectral_diversity/band_entropy.R`
- `Spectral_diversity/HSI_global_PCA.R`
- `Spectral_diversity/SA_entropy_bootstrapping.R`
- `Spectral_diversity/spectral_angle_entropy.R`
- `Spectral_diversity/spectral_angle_entropy_movingwindow.R`
- `Spectral_diversity/spectral_angle_entropy_raster.R`

### `scripts/3_Analysis/`

- `Analysis_1.R`
- `Analysis_PDF.R`
- `Bootstrap_Analysis.R`
- `global_vs_local.R`
- `linearity_test.R`
- `LLM.R`
- `modeling.R`
- `modeling_test.R`
- `ouliers.R`
- `random Forest.R`
- `scratch paper.R`
- `Spatial_random_Forest.R`
- `within_between_variation.R`
- `Currently not relevant/Topo_Analysis.R`

### `scripts/auxilary/`

- `kml_shp_conversion.R`
- `missing tiles.R`
- `RGB_Conversion.R`
- `sp_summary.R`
- `Stats.R`
- `vegindex summerization effects.R`

### `scripts/visuals/`

- `edge_quads.R`
- `interactrive_crown_radius.R`
- `shadow_mask.R`
- `smooth_spectra_compare.R`
- `Spatial Correlation.R`
- `spatial_corr_plot.R`
- `Spectral_Signatures.R`

### Remaining Review Candidate Scripts

- `scripts/3_Analysis/scratch paper.R`
- `scripts/3_Analysis/modeling_test.R`

## File-Type Summary

| Extension | Count | Notes |
|---|---:|---|
| `.xml` | 12,509 | Raster auxiliary metadata, mainly spectral products |
| `.hdr` | 12,419 | ENVI raster headers |
| no extension | 12,394 | ENVI raster data files |
| `.tif` | 1,163 | Raster imagery and derived products |
| `.R` | 66 | R workflow scripts, including root-level and ad hoc test scripts |
| `.shp/.dbf/.prj/.shx` | 42 each | Shapefile component sets |
| `.csv` | 10 | Tabular tree, taxonomy, spectral, biodiversity, and environmental data |
| `.adf` | 7 | ArcGIS raster/grid files |
| `.png` | 5 | Exported figures |
| `.docx` | 5 | Manuscript/workshop documents |
| `.pdf` | 13 | Reports and phylogeny output |
| `.Rhistory` | 4 before cleanup, 0 after cleanup | R session-history files removed |

## Cleanup Notes

- Removed `.Rhistory` session files on 2026-06-15.
- Removed authorized `old/`, `Outdated/`, and `Currently not relevant/` directories on 2026-06-15.
- The cleanup removed approximately 107.8 GB from `Quad_Spectra/old/`.
- Current spectral heterogeneity products are in `Indices_SHPs/Spectral_diversitySHPs/`.
- Spectra cropped to quadrat extents are in `Quad_Spectra/`.
