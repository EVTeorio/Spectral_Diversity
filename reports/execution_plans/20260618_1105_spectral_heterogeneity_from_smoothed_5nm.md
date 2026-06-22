# ExecPlan: Spectral Heterogeneity from Smoothed 5 nm Quadrat Spectra

Date: 2026-06-18

## Objective

Create per-quadrat spectral heterogeneity values for each analysis resolution from the current smoothed and 5 nm resampled quadrat spectra.

## Requested Task

- Pick up from `Quad_Spectra/*_smooth_5nm`.
- Use the current logic in `SA_entropy_bootstrapping.R`.
- Consider whether all pixels in each quadrat can reasonably be used.
- If all-pixel entropy is not computationally reasonable, fall back to the mean of bootstrap values.
- Preserve and document the shadow masking threshold.

## Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `RESEARCH_OBJECTIVES.md`
- `reports/project_state.md`
- `reports/directory_map.md`
- `reports/data_dictionary.md`
- `scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R`
- `scripts/2_Indices Creation/Spectral_diversity/spectral_angle_entropy.R`
- `scripts/1_Data Processing/smoothing.R`
- `scripts/1_Data Processing/resampling.R`

## Proposed Changes

- Refactor `SA_entropy_bootstrapping.R` to loop over 10 m, 20 m, and 50 m current smoothed 5 nm input folders.
- Add an exact all-pixel entropy path only when the number of pairwise spectral angles is below a configurable limit.
- Use a blockwise histogram calculation so entropy can be computed without repeatedly growing large angle vectors.
- Fall back to bootstrap replicates for larger quadrats and use `boot_mean` as the primary heterogeneity value.
- Keep shadow masking parameters explicit: band `563 nm`, threshold `0.0305476`, and sunlit pixels retained as values greater than the threshold under the current direction setting.
- Write per-scale CSV and shapefile outputs with fields for method, pixel counts, exact entropy when available, and bootstrap summaries when used.

## Validation Plan

- Run a lightweight `testthat` suite on helper functions using small matrices.
- Parse the updated script with R 4.2.3.
- Avoid running the full raster workflow in validation unless explicitly requested because the full per-scale calculation is potentially long-running.
- Update task, validation, and project-state Markdown notes.

## Risks

- Exact all-pixel entropy is quadratic in pixel count, so it is only reasonable for small masked quadrats.
- Full bootstrap processing across all 10 m, 20 m, and 50 m rasters may take substantial time.
- Existing 20 m historical products used older folder names, so downstream joins may need updated column names after the current workflow is run.
