# Multiscale Spectral-Biodiversity Analysis

Last updated: 2026-06-24

## Task

Created a reproducible downstream analysis and PDF reporting workflow for the current combined quadrat analysis tables.

## Inputs

- `quadrat_analysis_10m.csv`
- `quadrat_analysis_20m.csv`
- `quadrat_analysis_50m.csv`
- `scripts/3_Analysis/LLM.R` and `scripts/3_Analysis/Analysis_PDF.R` were used as style and method references.
- `RESEARCH_OBJECTIVES.md` supplied the scientific questions and hypothesis framing.

## Outputs

- `C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/Documents/PDFs/spectral_biodiversity_multiscale_findings.pdf`
- `C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/Documents/PDFs/spectral_biodiversity_model_appendix.pdf`
- `reports/figures/multiscale_spectral_biodiversity/`
- `reports/tables/multiscale_spectral_biodiversity/`

## Methods

- Primary response: `spec_sa`, the SA entropy mean. For most quadrats it is the mean of 70 bootstrap iterations using up to 5,000 retained pixels; the small exact subset uses all retained pixel pairs.
- Secondary responses: `spec_pca_mean`, `spec_rao_q`, and `spec_alpha`.
- Primary biodiversity predictors: `phy_faith`, `phy_afaith`, and `sp_shannon`.
- Environmental controls: `env_elev` and `env_tri11`.
- Predictors and responses were standardized within model datasets.
- Candidate OLS models were compared using adjusted R2 and AIC.
- Residual spatial autocorrelation was evaluated with 8-nearest-neighbor Moran's I.
- Documented 10 m and 20 m edge quadrats were excluded from primary analysis; 50 m used all quadrats because no separate 50 m edge rule is documented.

## Coverage Summary

| Scale | Total quadrats | Primary quadrats | Complete SA mean | Complete PCA mean | Edge flagged |
|---|---:|---:|---:|---:|---:|
| 10m | 2000 | 1656 | 1587 | 1440 | 344 |
| 20m |  500 |  414 |  405 |  360 |  86 |
| 50m |   80 |   80 |   80 |   74 |   0 |
