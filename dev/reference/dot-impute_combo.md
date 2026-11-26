# Mixed combo imputation: configurable MAR + MNAR methods

Two-stage imputation:

1.  MAR/MCAR: dispatched to mar_method

2.  MNAR: dispatched to mnar_method

## Usage

``` r
.impute_combo(
  x,
  condition,
  prop_na_in_condition = 1,
  prop_present_in_other_condition = 0,
  min_present_in_other_condition = 1,
  require_n_other_conditions = 1,
  mar_method = "Impseqrob",
  mnar_method = "min",
  method_args = list(),
  with_value = NA_real_
)
```

## Arguments

- x:

  Log2 intensity matrix

- condition:

  Condition vector aligned with columns

- prop_na_in_condition:

  Proportion of NA to classify as MNAR (default: 1.0)

- prop_present_in_other_condition:

  Proportion present in other conditions (default: 0.0)

- min_present_in_other_condition:

  Minimum present values (default: 1)

- require_n_other_conditions:

  Required conditions with presence (default: 1)

- mar_method:

  MAR imputation method (default: "Impseqrob")

- mnar_method:

  MNAR imputation method (default: "min")

- method_args:

  Named list of per-method argument lists

- with_value:

  Constant for method "with"

## Value

List with x_imputed, mnar_mask, mar_mask, summary
