# Project Governance

## Working Directory

`C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity`

## Operating Rules

- Use R 4.2.3-compatible code. Review nearby scripts before changing R workflows so naming, loops, and output patterns remain consistent.
- Update `reports/directory_map.md`, `reports/data_dictionary.md`, and `reports/project_state.md` whenever a material directory, dataset, or workflow decision changes.
- Keep current derived data in `Quad_Values/`, analytical scripts in `scripts/`, and dated manuscript work in `Documents/Paper/`.
- Never modify `Documents/`, manuscript-state files, or PDFs unless the user explicitly requests it.
- Never delete, move, or overwrite files outside the stated task without approval. For cleanup, first identify exact paths, establish whether files are duplicates or archival variants, and record the decision.
- Do not commit, push, rebase, alter system paths, or alter system environment variables.

## Current Data and Scale Rules

- Confirmed base quadrat spectra are `Quad_Spectra/10m/`, `20m/`, and `50m/`.
- `Quad_Spectra/10m_smooth_5nm/` is the intended current 10 m resampled input but was accidentally deleted and will be restored by the user.
- Intermediate `*_smooth/` directories were deliberately removed to conserve space and may be regenerated when needed.
- `Quad_Spectra/25m/`, `25m_smooth_5nm/`, and `33m/` support planned scale expansion. Empty or partial folders at these scales are intentional and must not be removed as cleanup artifacts.
- Current primary spectral metrics are original-PCA and vector-normalized PCA mean distance plus 2D alpha-hull area. Current biodiversity metrics are the species and phylogenetic summaries in `Quad_Values/Diversity_SHPs/` and the combined analysis tables.
- Preserve both PCA `.rds` basis objects and their variance-explained tables. Convex-hull, 3D hull, spectral-angle entropy, bootstrap, and sensitivity products are retained as archive/provenance data unless explicitly reopened. They are not current primary metrics.
- `Quad_Values/Spectral_diversitySHPs/` is the canonical home for current combined spectral-heterogeneity summaries and PCA basis objects. Root-level copies of those combined summary CSVs are redundant copies.

## Documentation and Planning

For complex changes, create a dated execution plan in `reports/execution_plans/` before editing. The plan must state objective, scope, files reviewed, proposed changes, validation, and risks.

Maintain documentation as a concise current-state record:

- `reports/project_state.md`: current objective, recent direction changes, active inputs, primary outputs, pending decisions, and known issues.
- `reports/data_dictionary.md`: canonical datasets, locations, schemas/keys, status, and archived or duplicate variants.
- `reports/directory_map.md`: live directory structure and intended role of each major location.

Historical task, validation, analysis, and execution-plan Markdown files are provenance. Do not delete them during routine documentation consolidation unless explicitly authorized.

## Cleanup Rules

- In `reports/`, remove only figures or tables that are verified duplicates. Keep the latest suitable rendering when an iteration is demonstrably superseded.
- Do not delete a similarly named plot or table merely because it covers related metrics; verify byte identity, source data, report references, and analysis purpose first.
- Retain compact archival tables when they preserve reproducibility or explain a prior analytical decision. Prefer one canonical copy when files are byte-identical and no scripts or current reports require the duplicate path.
- Existing folders named `old`, `Outdated`, and `Currently not relevant` were previously authorized cleanup targets; do not extend that authorization to new folders without review.

## Validation

For code or data changes, run focused validation proportional to risk and document the result. For documentation-only work, verify links/paths, file inventories, protected locations, and any stated duplicate decisions.

Before final completion of a substantial task, run an interactive R completion beep with `beepr` when R access is available. Report if `beepr` is unavailable.
