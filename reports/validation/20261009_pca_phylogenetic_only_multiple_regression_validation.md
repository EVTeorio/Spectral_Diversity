# PCA Phylogenetic-Only Multiple Regression Validation

Date: 2026-10-09

## Validation Checks

| check | expected | observed | passed |
| --- | --- | --- | --- |
| full_model_row_count | 18 | 18 | TRUE |
| coefficient_row_count | 72 | 72 | TRUE |
| partial_r2_row_count | 54 | 54 | TRUE |
| vif_row_count | 54 | 54 | TRUE |
| all_full_models_estimable | 18 | 18 | TRUE |
| finite_full_model_r2 | 18 | 18 | TRUE |
| complete_case_n_positive | 18 | 18 | TRUE |
| paper_table_created | exists |  1 | TRUE |

## Notes

- All models use complete cases for the response and the three phylogenetic predictors.
- Complete-case counts match PCA metric availability from the combined quadrat tables.
