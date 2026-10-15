# Directory Map

Last updated: 2026-10-14

```text
Spectral_Diversity/
    CODEX_AGENT_GUIDELINES.md       governance and cleanup rules
    RESEARCH_OBJECTIVES.md          current scientific and analytical focus
    quadrat_analysis_*m.csv         joined analysis-ready tables for established scales

    Documents/                      protected manuscript, paper, and PDF materials; do not modify in routine cleanup
    Functions/                      legacy and support R functions
    HSI_NA_trimmed/                 trimmed source hyperspectral imagery
    Quad_Scale_SHPs/                10 m, 20 m, 25 m, 33 m, and 50 m quadrat boundaries
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
| `Quad_Spectra/10m_smooth_5nm` | 10 m current resampled input | Accidentally deleted; user will restore. |
| `Quad_Spectra/20m_smooth_5nm`, `50m_smooth_5nm` | Current resampled inputs | Current. |
| `Quad_Spectra/25m`, `25m_smooth_5nm`, `33m` | Scale-expansion resources | Planned/partial; retain. |
| `Quad_Values/Diversity_SHPs/` | Species and phylogenetic diversity values | Canonical. |
| `Quad_Values/Enviro_SHPs/` | Environmental covariates | Canonical. |
| `Quad_Values/Spectral_diversitySHPs/` | Current combined spectral summaries and PCA bases | Canonical. |
| `reports/figures/`, `reports/tables/` | Derived analysis artifacts | Remove only verified duplicates; retain provenance otherwise. |

## Cleanup Boundaries

- `Documents/`, manuscript-state materials, and PDFs are out of scope for routine cleanup.
- Intermediate `*_smooth/` spectral directories were intentionally removed as recreatable products.
- Root-level combined spectral summary CSVs duplicate the canonical copies in `Quad_Values/Spectral_diversitySHPs/`; do not create new root-level copies.
- Do not remove empty or partial 25 m and 33 m locations; they are planned analysis scales.
