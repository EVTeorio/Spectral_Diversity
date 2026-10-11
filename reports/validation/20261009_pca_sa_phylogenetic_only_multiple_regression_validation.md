# PCA and Spectral Angle Phylogenetic-Only Multiple Regression Validation

Date: 2026-10-09

## Validation Checks

| check | expected | observed | passed |
| --- | --- | --- | --- |
| full_model_row_count | 21 | 21 | TRUE |
| coefficient_row_count | 84 | 84 | TRUE |
| partial_r2_row_count | 63 | 63 | TRUE |
| vif_row_count | 63 | 63 | TRUE |
| all_full_models_estimable | 21 | 21 | TRUE |
| finite_full_model_r2 | 21 | 21 | TRUE |
| complete_case_n_positive | 21 | 21 | TRUE |
| paper_table_created | exists |  1 | TRUE |

## Notes

- All models use complete cases for the response and the three phylogenetic predictors.
- Complete-case counts match PCA metric availability from the combined quadrat tables.
