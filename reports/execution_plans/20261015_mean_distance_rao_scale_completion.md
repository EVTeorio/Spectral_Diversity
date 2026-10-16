# Mean-distance and Rao's Q completion and restored 10 m verification

Date: 2026-10-15

## Objective and scope

Verify the restored 10 m smoothed 5 nm spectra by recomputing vector-normalized
PC1-PC3 mean distance and comparing every available eligible quadrat with saved
values. Fill the 25 m original and vector-normalized PCA mean-distance and Rao's Q
outputs. Audit 10 m, 20 m, 25 m, and 50 m coverage. Do not calculate alpha hulls or
process 33 m or 100 m. Do not modify manuscript assets or saved PCA bases.

## Reviewed inputs

- CODEX_AGENT_GUIDELINES.md and current project-state, dictionary, and map.
- scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R
- scripts/3_Analysis/alpha_hull_all_pixels_values.R (projection and mask reference).
- Existing canonical spectral summaries, component tables, and PCA basis objects.
- Restored 10m_smooth_5nm and prepared 25m_smooth_5nm ENVI inputs and scale boundaries.

## Proposed changes

1. Add a focused R workflow that reuses current masking, normalization, and saved
   PCA projection functions without invoking hull calculations or rebuilding PCA.
2. Run preflight checks on IDs, raster bands, saved basis metadata, and known
   metric coverage. Compare a pilot restored raster before the complete rerun.
3. Save independent 10 m verification results and per-quadrat differences; retain
   the historical values as the reference rather than overwriting them.
4. Generate 25 m distance/Rao component tables and a canonical distance/Rao-only
   summary with full-boundary coverage and explicit missing-raster statuses.
5. Per the user's subsequent clarification, calculate all available 25 m rasters
   without adding exclusion flags; the user will handle exclusions manually.
6. Save coverage and validation reports and update the three core state documents.

## Validation

- Confirm PC1-PC3 and all retained pixels, the saved 563 nm shadow threshold,
  vector normalization, and squared-distance Rao formula.
- Compare retained-pixel counts and finite values by ID, with explicit absolute
  and relative tolerances, duplicate/missing IDs, and maximum differences.
- Verify existing component and canonical mean/Rao fields agree at established scales.
- Validate output IDs against boundaries and distinguish missing rasters,
  exclusions, processing failures, and finite metrics.
- Verify PCA basis and existing summary hashes remain unchanged.
- Run the interactive completion beep when the installed beepr library is accessible.

## Risks and boundaries

- Restored files may differ in content or may be incomplete; report discrepancies
  without replacing the saved reference.
- 25 m exclusion decisions belong to the user; no new exclusion flags or spatial
  exclusion calculations are part of this task.
- Missing imagery is not repaired by inventing values or regenerating other scales.
- Existing user R packages require execution with access to the user library.
