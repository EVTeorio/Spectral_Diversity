# Task Report: Bootstrap Variation Analysis

Date: 2026-06-18

## Objective

Compare within-quadrat bootstrap variation against between-quadrat spectral heterogeneity variation to assess whether the current bootstrap-derived spectral entropy values are stable enough for downstream analysis.

## Inputs

- `Indices_SHPs/10m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/20m_SA_entropy_smooth_masked_5nm_summary.csv`
- `Indices_SHPs/50m_SA_entropy_smooth_masked_5nm_summary.csv`
- matching `_boot_wide.csv` files for 10 m, 20 m, and 50 m

## Script Added

- `scripts/3_Analysis/bootstrap_variation_analysis.R`

## Outputs Created

- `reports/analysis/20260618_bootstrap_variation_analysis.md`
- `reports/tables/bootstrap_variation/bootstrap_variation_scale_summary.csv`
- `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv`
- `reports/tables/bootstrap_variation/bootstrap_variation_top_unstable_quadrats.csv`
- `reports/figures/bootstrap_variation/within_vs_between_sd.png`
- `reports/figures/bootstrap_variation/bootstrap_cv_distribution.png`
- `reports/figures/bootstrap_variation/bootstrap_sd_vs_entropy.png`
- `reports/figures/bootstrap_variation/spectral_entropy_histograms_by_scale.png`
- `reports/figures/bootstrap_variation/pixel_count_vs_bootstrap_cv.png`
- `reports/figures/bootstrap_variation/variance_partition_diagnostic.png`
- `reports/figures/bootstrap_variation/top_unstable_quadrat_bootstrap_distributions.png`

## Key Results

- The bootstrap means are usable as primary spectral heterogeneity estimates when quality-control fields are carried forward.
- Raw within-quadrat bootstrap variation is not uniformly small enough to ignore.
- Median bootstrap CV is low at all scales:
  - 10 m: 0.4%
  - 20 m: 0.5%
  - 50 m: 0.5%
- Some quadrats show elevated bootstrap CV and should be flagged:
  - 10 m: 372 of 1,866 bootstrap-estimated quadrats exceed 5% CV; 127 exceed 10% CV.
  - 20 m: 60 of 477 exceed 5% CV; 12 exceed 10% CV.
  - 50 m: 11 of 79 exceed 5% CV; 0 exceed 10% CV.
- The 95th percentile bootstrap CV across all bootstrap-estimated quadrats is 10.4%.
- The spectral entropy histogram figure shows mild right skew at 10 m, near-symmetry at 20 m, and stronger right skew at 50 m.

## Recommendation

Use `spectral_entropy` / `boot_mean` as the main heterogeneity estimate, but carry `boot_sd`, `boot_cv`, and `method` into downstream tables. Run sensitivity checks that flag or exclude high-CV quadrats, especially above 5% or 10% bootstrap CV.
