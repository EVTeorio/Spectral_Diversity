# Directory Map

Last updated: 2026-10-15

```text
Spectral_Diversity/
    CODEX_AGENT_GUIDELINES.md       governance and cleanup rules
    RESEARCH_OBJECTIVES.md          current scientific and analytical focus
    quadrat_analysis_*m.csv         joined analysis-ready tables for established scales

    Documents/                      protected manuscript, paper, and PDF materials; do not modify in routine cleanup
    Functions/                      legacy and support R functions
    HSI_NA_trimmed/                 trimmed source hyperspectral imagery
    Quad_Scale_SHPs/                10 m, 20 m, 25 m, 33 m, 50 m, and 100 m quadrat boundaries
    Quad_Spectra/                   base spectra, resampled spectra, planned scale folders, and visual derivatives
    Quad_Values/                    canonical derived spectral, diversity, and environmental values
    reports/                        current-state documentation plus analysis, figure, table, task, plan, and validation provenance
    scripts/                        data processing, index creation, analysis, utility, and visual R scripts
    tests/                          formal tests and lightweight reproduced-raster validation artifacts
    tools/                          project support utilities
    zotero_index/                   local literature-index support files
    logs/                           processing logs
    Visuals/                        exported visual products
```

## Core Locations

| Location | Role | Current status |
|---|---|---|
| `Documents/` | Manuscript and PDF materials | Protected from this consolidation. |
| `Quad_Spectra/10m`, `20m`, `50m` | Confirmed base spectral inputs | Current. |
| `Quad_Spectra/10m_smooth_5nm` | 10 m current resampled input | Restored from the data repository; 1,909 raster IDs match the saved spectral footprint. |
| `Quad_Spectra/20m_smooth_5nm`, `50m_smooth_5nm` | Current resampled inputs | Current. |
| `Quad_Spectra/25m`, `25m_smooth_5nm` | Active scale expansion | 312 prepared rasters; 320 boundary cells. |
| `Quad_Spectra/33m`, `Quad_Scale_SHPs/PR_33m.*`, `PR_100m.*` | Deferred scale resources | No 33 m or 100 m processing in the current task; retain. |
| `Quad_Values/Diversity_SHPs/` | Species and phylogenetic diversity values | Canonical. |
| `Quad_Values/Enviro_SHPs/` | Environmental covariates | Canonical. |
| `Quad_Values/Spectral_diversitySHPs/` | Current combined spectral summaries and PCA bases | Canonical. |
| `Quad_Values/Spectral_diversitySHPs/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv` | 25 m distance/Rao summary in both bases | 320 grid rows; 312 finite values per metric/basis; no hulls or new exclusions. |
| `reports/tables/20261015_pca_mean_rao_completion/` | Restored 10 m verification and 25 m completion audit | Per-quadrat comparisons, coverage, validation, and resumable computation chunks. |
| `reports/tables/20261015_alpha_half_small_scales/` | Authorized alpha = 0.5 hull run | Per-quadrat checkpoints, timing, validation, and status summaries for 10 m, 20 m, and 25 m only. |
| `Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_*m_all_pixels.csv` | Canonical all-pixel alpha = 0.5 areas | Completed 10 m, 20 m, and 25 m tables in both bases, with full-grid statuses. Only 10 m original-PCA 1000_c remains unresolved. User-authorized manual orientation repairs: 20 m 413; 10 m 1402_b, 1019_a, and 1816_c. |
| `reports/tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/` | Manual repair provenance for original-PCA 413 | Original-table backups, signed arc traversal, area audit, and focused validation; original saved hull remains in the adjacent retry directory. |
| `reports/tables/20261015_alpha_half_small_scales/manual_orientation_three_10m/` | Reviewed manual repairs of 1402_b, 1019_a, and 1816_c | Before-table backups, signed traversals, audit, hashes and validation. Saved original hulls remain in `diagnose_four_10m/`. |
| `reports/figures/`, `reports/tables/` | Derived analysis artifacts | Remove only verified duplicates; retain provenance otherwise. |

## Cleanup Boundaries

- `Documents/`, manuscript-state materials, and PDFs are out of scope for routine cleanup.
- Intermediate `*_smooth/` spectral directories were intentionally removed as recreatable products.
- Root-level combined spectral summary CSVs duplicate the canonical copies in `Quad_Values/Spectral_diversitySHPs/`; do not create new root-level copies.
- Do not remove empty or partial 25 m and 33 m locations; they are planned analysis scales.
