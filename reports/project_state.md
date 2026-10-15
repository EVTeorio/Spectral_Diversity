# Project State

Last updated: 2026-10-14

## Current Objective

Maintain a compact, reproducible spectral-biodiversity workflow while expanding beyond 10 m, 20 m, and 50 m to planned 25 m and 33 m scales.

## Active Analysis

- Primary spectral metrics: original-PCA and vector-normalized PCA mean distance and 2D alpha-hull area.
- Biodiversity metrics: richness, Shannon, Simpson, evenness, Faith PD, phylogenetic Rao's Q, and abundance-weighted Faith PD.
- Canonical spectral summaries: `Quad_Values/Spectral_diversitySHPs/`.
- Canonical biodiversity summaries: `Quad_Values/Diversity_SHPs/`.
- Current joined inputs: `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.

## Recent Direction Changes

- The workflow now prioritizes both original-PCA and vector-normalized PCA distance and alpha-hull area over spectral angle entropy, convex hulls, and 3D hull products.
- Historical metric products are retained as compact provenance, not deleted by default.
- Intermediate `*_smooth/` spectra were intentionally removed to conserve space because they can be regenerated.
- `10m_smooth_5nm` was accidentally deleted and will be restored by the user.
- 25 m and 33 m are planned scales; their incomplete resources are intentional.
- Repository consolidation must not alter `Documents/`, manuscript-state assets, or PDFs.

## Completed Work

- Established current diversity, environmental, and spectral output locations.
- Generated and documented 10 m, 20 m, and 50 m biodiversity and environmental summaries.
- Produced global original and vector-normalized PCA bases and spectral-summary tables.
- Created joined analysis tables for the established scales.
- Completed historical sensitivity, bootstrap, correlation, regression, and alpha-hull exploratory work retained under `reports/`.
- Refreshed governance, research objectives, directory mapping, and the data dictionary to distinguish current products from archived provenance.

## Pending Work

- Restore `Quad_Spectra/10m_smooth_5nm/` before rerunning any 10 m spectral workflow.
- Build and validate 25 m and 33 m scale spectra, diversity, environmental outputs, and joined tables.
- Decide whether the three root-level duplicated combined spectral-summary CSVs should be removed after downstream reference checks.
- Decide which compact archive products are worth retaining once 25 m and 33 m outputs exist.
- Apply spatial/environmental sensitivity work before final inferential claims where needed.

## Known Issues and Decisions Needed

- `Quad_Values/` contains historical 20 m-only outputs that are not current inputs but may be useful provenance. The earlier `Spectral_diversitySHPs_first/` archive has been removed; current canonical spectral products remain in `Spectral_diversitySHPs/`.
- The root-level and `Spectral_diversitySHPs/` combined spectral-summary CSVs are exact duplicates; the latter is canonical.
- Preserve both PCA `.rds` basis objects and their variance-explained tables. Future modeling should explicitly select `spec_pca_mean`, `spec_alpha`, `spec_spca_mean`, and `spec_spca_alpha`.
- No byte-identical figures were found in `reports/figures/`; do not remove figures based solely on topic or filename similarity.
