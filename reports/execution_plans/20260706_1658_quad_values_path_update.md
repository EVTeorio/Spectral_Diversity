# ExecPlan: Quad_Values Path Update

## Objective

Make `Quad_Values/` the canonical derived-output directory throughout scripts, reports, tables, and governance documentation.

## Requested Task

- Replace remaining `Quad_Values` references with `Quad_Values`.
- Update R scripts that read or write derived spectral, diversity, environmental, and topographic outputs.
- Update Markdown and CSV artifacts so the documented paths match the renamed directory.
- Regenerate/update user-facing documentation where practical.

## Files To Review

- `reports/project_state.md`
- `reports/directory_map.md`
- `reports/data_dictionary.md`
- `CODEX_AGENT_GUIDELINES.md`
- R scripts under `scripts/`
- Markdown reports under `reports/`
- Sample-size design CSVs under `reports/tables/sample_size_effects/`

## Proposed Changes

- Perform a repository-wide text replacement from `Quad_Values` to `Quad_Values` across text artifacts.
- Manually revise narrative notes that previously framed `Quad_Values/` as a mismatch or pending migration.
- Validate that no text references to `Quad_Values` remain.
- Parse key R scripts after the rename where feasible.

## Validation Plan

- Search the repository for any remaining `Quad_Values` references.
- Confirm key output directories exist under `Quad_Values/`.
- Run lightweight R parse checks for modified R scripts.
- Run the completion beep with `beepr`.

## Risks

- Historical reports will now reflect the renamed current path rather than the original path used at the time those reports were created.
- Some generated binary artifacts, especially `.docx`, may require regeneration if they embed the old path.
