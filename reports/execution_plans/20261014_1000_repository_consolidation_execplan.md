# Repository Consolidation and Documentation Refresh

Date: 2026-10-14

## Objective

Refresh the repository's core operational documentation, identify duplicated or scattered derived objects, and remove only verified duplicate figures and tables under `reports/` while preserving current manuscript materials and recoverable archival data.

## Scope

- Update `CODEX_AGENT_GUIDELINES.md`, `RESEARCH_OBJECTIVES.md`, `reports/data_dictionary.md`, `reports/directory_map.md`, and `reports/project_state.md`.
- Audit duplicate or overlapping datasets in `Quad_Values/`, root-level analysis tables, and `reports/tables/`.
- Compare `reports/figures/` and `reports/tables/` by content, schema, provenance, and recency before removing only confirmed duplicates.
- Preserve `Documents/`, all manuscript-state materials, and PDFs without modification.

## Current Decisions

- `Quad_Spectra/10m_smooth_5nm/` was accidentally deleted and will be restored by the user; it remains the intended current 10 m analysis input.
- Intermediate `*_smooth/` directories were deliberately removed to save space and may be regenerated.
- `25m` and `33m` are planned analytical scales; empty or partial folders at those scales are intentional placeholders, not cleanup targets.
- Current priority spectral metrics are vector-normalized PCA mean distance and 2D alpha-hull area. Current biodiversity metrics are the current species and phylogenetic diversity summaries.
- Older spectral metrics and intermediate products are retained only as compact archive/provenance data unless they are exact duplicates.

## Files To Review

- Core governance and research Markdown documents.
- `scripts/3_Analysis/sample_2d_pca_mean_rao_values.R`
- `scripts/3_Analysis/alpha_hull_all_pixels_values.R`
- `scripts/3_Analysis/combine_quadrat_analysis_tables.R`
- Current CSVs and shapefiles under `Quad_Values/`.
- Figures and tables under `reports/`.

## Proposed Changes

- Condense core project-state and research documentation around current inputs, current metrics, current scales, historical direction changes, and archive rules.
- Update the data dictionary and directory map to reflect the live tree and the deliberately absent/recoverable intermediates.
- Write a compact consolidation inventory identifying canonical datasets, duplicates, archived products, and candidates requiring user direction.
- Delete only report figures/tables that are content-identical or demonstrably superseded duplicate renderings; retain the newest suitable version.

## Validation Plan

- Verify protected `Documents/` and PDF files are unchanged.
- Record pre/post counts and paths for any removed report artifacts.
- Verify documentation references the current 10 m, 20 m, 25 m, 33 m, and 50 m scale status accurately.
- Confirm all deletion candidates were located under `reports/figures/` or `reports/tables/` and met the duplicate criteria.

## Risks

- Semantically similar plots can communicate different analyses; filename similarity alone is insufficient for deletion.
- Existing tables contain both current values and provenance/diagnostic variants; compact archival preservation is preferable to broad deletion.
- Some live folders are intentionally incomplete during scale expansion and must not be classified as erroneous.
