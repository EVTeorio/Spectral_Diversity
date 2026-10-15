# Task Report: Repository Consolidation

Date: 2026-10-14

## Completed

- Refreshed the governance, research objective, data dictionary, directory map, and project-state documents around the current PCA distance/area focus.
- Documented the planned 25 m and 33 m expansion, deliberate removal of recreatable `*_smooth/` intermediates, and pending restoration of `10m_smooth_5nm/`.
- Classified current, archived, and redundant data objects without modifying protected `Documents/` or PDF materials.
- Removed one exact duplicate report table; its canonical copy remains at `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv`.

## Canonical Decision

The canonical combined spectral summaries are the three `spectral_heterogeneity_*m_smooth_masked_5nm_summary.csv` files in `Quad_Values/Spectral_diversitySHPs/`. Byte-identical root-level copies remain in place pending an explicit data-deletion decision.

## Protected Scope

No files under `Documents/` and no PDF files were changed.
