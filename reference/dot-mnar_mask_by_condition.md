# MNAR mask by condition with evidence in other conditions

Identifies MNAR candidate cells: NA in one condition but with sufficient
presence in other conditions.

## Usage

``` r
.mnar_mask_by_condition(
  x,
  condition,
  prop_na_in_condition = 1,
  prop_present_in_other_condition = 0,
  min_present_in_other_condition = 1,
  require_n_other_conditions = 1,
  drop_empty_levels = TRUE
)
```

## Arguments

- x:

  Intensity matrix (proteins x samples)

- condition:

  Condition vector aligned with columns

- prop_na_in_condition:

  Minimum proportion of NA in the condition (default: 1.0 = 100%)

- prop_present_in_other_condition:

  Minimum proportion of non-NA in other condition (default: 0.0)

- min_present_in_other_condition:

  Minimum absolute non-NA count in other condition (default: 1)

- require_n_other_conditions:

  Number of other conditions that must meet criteria (default: 1)

- drop_empty_levels:

  Drop empty factor levels (default: TRUE)

## Value

Logical matrix of same dimensions indicating MNAR cells
