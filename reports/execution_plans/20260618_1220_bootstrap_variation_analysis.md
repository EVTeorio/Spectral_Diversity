# ExecPlan: Bootstrap Variation Analysis

Date: 2026-06-18

## Objective

Analyze the spectral angle entropy bootstrap outputs to compare within-quadrat bootstrap variation against between-quadrat spectral heterogeneity variation.

## Requested Task

Create an analysis report with graphs and visuals to determine whether within-quadrat bootstrap variation is small enough for the spectral heterogeneity estimates to be considered valid.

## Inputs

- `Indices_SHPs/10m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/20m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/50m_SA_entropy_smooth_masked_5nm_summary.csv`
- matching `_boot_wide.csv` files for each scale

## Proposed Changes

- Add an R analysis script under `scripts/3_Analysis/`.
- Generate tabular summaries comparing:
  - bootstrap SD and CV within quadrats
  - between-quadrat SD and variance of spectral entropy
  - ratios of within to between variation
  - standard error of bootstrap means
  - outlier and instability diagnostics
- Generate figures under `reports/figures/bootstrap_variation/`.
- Write a Markdown analysis report under `reports/analysis/`.

## Validation Plan

- Confirm all expected input files exist and have 70 bootstrap columns.
- Confirm generated figures exist and are non-empty.
- Confirm summary tables have one row per scale.
- Review key results for missing or extreme values.

## Risks

- Some quadrats have insufficient pixels or missing bootstrap replicates and should be excluded from within-bootstrap summaries.
- Bootstrap SD measures variation among subsampled estimates, while the uncertainty of the reported `boot_mean` is better represented by `boot_sd / sqrt(n_boot)`.
- The bootstrap implementation uses sampled pixel pairs for large pixel subsamples, so the report should interpret results as computational bootstrap stability rather than a complete census of all possible pixel pairs.
