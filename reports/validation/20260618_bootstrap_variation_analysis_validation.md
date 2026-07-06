# Validation Report: Bootstrap Variation Analysis

Date: 2026-06-18

## Validation Scope

Validate the bootstrap variation analysis report, figures, and diagnostic tables generated from current spectral heterogeneity bootstrap outputs.

## Checks Performed

| Check | Result | Notes |
|---|---|---|
| Input files exist | Pass | Summary and wide bootstrap CSVs were read for 10 m, 20 m, and 50 m. |
| Bootstrap columns | Pass | Each wide bootstrap file contains 70 bootstrap replicate columns. |
| Analysis script execution | Pass | `scripts/3_Analysis/bootstrap_variation_analysis.R` completed successfully. |
| Tables generated | Pass | Scale summary, quadrat diagnostics, and top-unstable-quadrat tables were written under `reports/tables/bootstrap_variation/`. |
| Figures generated | Pass | Seven PNG figures were written under `reports/figures/bootstrap_variation/`. |
| Figure readability | Pass | Main figures were visually inspected; PNGs were regenerated with white backgrounds for readability. The spectral entropy histogram was visually inspected after generation. |
| Report generated | Pass | `reports/analysis/20260618_bootstrap_variation_analysis.md` was written with tables, figure links, interpretation, and recommendation. |

## Main Finding

The bootstrap means are usable with quality-control flags, but raw within-quadrat bootstrap variation is not uniformly small enough to ignore.

## Quantitative Summary

| Scale | Median Bootstrap CV | P90 Bootstrap CV | Bootstrap Quads > 5% CV | Bootstrap Quads > 10% CV | Stability Class |
|---|---:|---:|---:|---:|---|
| 10 m | 0.4% | 8.5% | 372 of 1,866 | 127 of 1,866 | Review |
| 20 m | 0.5% | 5.7% | 60 of 477 | 12 of 477 | Review |
| 50 m | 0.5% | 5.8% | 11 of 79 | 0 of 79 | Review |

## Distribution Shape

The added spectral heterogeneity histogram indicates mild right skew at 10 m, near-symmetry at 20 m, and stronger right skew at 50 m. Numeric skewness values in the scale summary table are approximately 0.275 for 10 m, -0.002 for 20 m, and 1.187 for 50 m.

## Recommendation

Use `spectral_entropy` / `boot_mean` as the primary heterogeneity mean, but carry `boot_sd`, `boot_cv`, bootstrap-mean CI fields, and `method` forward into downstream analysis tables. Include sensitivity checks that flag or exclude high-CV or wide-CI quadrats.
