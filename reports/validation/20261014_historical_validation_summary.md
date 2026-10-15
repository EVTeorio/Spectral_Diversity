# Historical Validation Summary

Last updated: 2026-10-14

This summary replaces completed validation notes that supported earlier workflow stages. It preserves the critical validation conclusions needed to interpret the current data.

## Key Validations Retained in Summary

- Current plant-diversity and environmental tables were validated at 10 m, 20 m, and 50 m using scale plus quadrat ID joins.
- The original and vector-normalized PCA bases were validated from masked 10 m spectra. Both PCA `.rds` objects and variance-explained tables remain preserved.
- Combined spectral summaries and joined quadrat tables were checked for row counts, identifiers, missingness, and expected PCA fields.
- Spectral-angle entropy bootstrap outputs were validated as reproducible estimates with quality-control fields, but are not current primary response metrics.
- Alpha-hull runs record manual exclusions and non-finite geometry statuses explicitly rather than silently omitting them.
- Historical simple-regression and OLS results require caution because spatial autocorrelation was identified in several priority relationships.

## Current Detailed Validation Records

Detailed records remain for the current manuscript-direction and October PCA/alpha-hull work, including the October 9 regressions, October 10 benchmark, October 11 sampled PCA run, and October 14 repository consolidation.
