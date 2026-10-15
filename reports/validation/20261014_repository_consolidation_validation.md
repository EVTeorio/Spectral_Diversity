# Validation Report: Repository Consolidation

Date: 2026-10-14

| Check | Result | Notes |
|---|---|---|
| Protected manuscript scope | Passed | No files under `Documents/` or PDFs were changed. |
| 25 m and 33 m status | Passed | Core documents identify these as planned scales, not cleanup targets. |
| Missing 10 m resampled spectra | Passed | Core documents record that `Quad_Spectra/10m_smooth_5nm/` will be restored by the user. |
| Current metric focus | Passed | Core documents identify vector-normalized PCA distance and 2D alpha-hull area as priority metrics. |
| Current biodiversity focus | Passed | Core documents identify the current species and phylogenetic diversity metrics. |
| Spectral-summary duplication | Passed | Three root-level combined summaries were verified byte-identical to canonical copies in `Quad_Values/Spectral_diversitySHPs/`; no data copies were deleted. |
| Report figures | Passed | No byte-identical figures were found; none were removed. |
| Report tables | Passed | One unreferenced exact duplicate SA sample-size design table was removed; canonical copy remains at `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv`. |

## Remaining Decision

The three root-level duplicated combined spectral-summary CSVs should be removed only after the user confirms that the canonical `Quad_Values/Spectral_diversitySHPs/` paths should become the sole retained copies.
