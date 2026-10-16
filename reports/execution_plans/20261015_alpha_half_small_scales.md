# All-pixel alpha = 0.5 hulls at smaller scales

Date: 2026-10-15

## Objective and authorized scope

Complete original and vector-normalized PC1-PC2 alpha-hull areas at alpha = 0.5
at 10 m, then 20 m, then 25 m. Stop before 50 m. Do not process 33 m or 100 m.
Reuse validated finite existing alpha = 0.5 values and retry failed small-scale
values. Preserve both saved PCA bases, previous mean/Rao outputs, and historical
hull products. The user handles exclusions manually; do not add exclusion flags.
Keep existing 10 m/20 m reference exclusions for continuity, with no new ones.

## Reviewed files

- CODEX_AGENT_GUIDELINES.md and reports/project_state.md.
- scripts/3_Analysis/alpha_hull_all_pixels_values.R.
- scripts/3_Analysis/complete_pca_mean_rao_and_verify_10m.R.
- reports/analysis/20261012_alpha_hull_all_pixels_values_multiscale_alpha0_5_1_1_5.md.
- Existing alpha-hull results, canonical spectral summaries, saved PCA bases,
  and completed 25 m distance/Rao outputs.

## Changes

1. Add a focused resumable workflow restricted to 10 m, 20 m, and 25 m.
2. Reuse the established shadow mask, cleaning, normalization, and PCA projection.
3. Calculate R alphahull geometry from all unique retained PC1-PC2 coordinates;
   do not sample pixels, change alpha, manually jitter coordinates, or substitute
   a hull. Preserve warnings from the installed triangulator's native retries.
4. Save each result with provenance, runtime, pixel counts, and explicit status;
   write canonical per-scale alpha-only tables and an audit/status report.
5. Update current project-state, map, and dictionary after validation.

## Validated recovery for package area-evaluation failures

The installed package's endpoint-reordering area routine can fail with
`the condition has length > 1`. An exact Green's-theorem line integral of the
same oriented circular arcs was derived and validated: 100 synthetic cases
(clouds, disconnected components, holes, translations), three finite 10 m
comparisons, and four actual 20 m/25 m basis comparisons agree within 1e-8.
The largest synthetic difference is approximately 2.12e-12. A failed 10 m area
was recovered without changing geometry. Use this only after the main run for
failed evaluations, preserve first-attempt CSVs, and record provenance. Successful
package results remain unchanged. Hull-construction failures may remain unresolved.

An equivalent endpoint/chord integral was subsequently validated on 100 synthetic
cases (maximum discrepancy approximately 7.38e-14) and actual hulls. Exact vertex-ID
balance, circle consistency, and a convex-area bound were required. This recovered
two additional cases with closed stored endpoints. The endpoint form also
recovered the remaining 25 m original area following a native triangulation retry.
Final outcome: 13 recovered cases,
five unresolved original-PCA cases, all vector-normalized target areas finite, and
no 50 m processing. Audit records preserve repeated attempts and method provenance.

## Validation

- Reuse keys include scale, quadrat ID, PCA basis, and alpha = 0.5.
- Confirm unique keys and full intended input/expected-ID coverage.
- Check band alignment, finite nonnegative areas, and retained-pixel counts
  against the established distance/Rao pipeline.
- Confirm original reference CSVs, both saved bases, and all 50 m results remain
  unchanged by hashes. Preserve failures explicitly rather than mark them complete.
- Record completion notification via beepr when available.

## Risks

All-pixel geometry can be slow or numerically nonfinite. Longer area calculation
timeouts may recover previous timeouts, but unresolved cases will remain explicit.
Missing rasters cannot receive inferred areas. Existing alpha = 1 fields must not
be silently relabeled or replaced by the new alpha = 0.5 values.
