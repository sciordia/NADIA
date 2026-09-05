# Pre-filter proteins by MNAR rules (for combo mode)

Keeps a protein if: (A) It has MNAR in \>= 1 condition, OR (B) It has
sufficient signal in \>= require_n_other_conditions conditions

## Usage

``` r
.prefilter_by_rules(
  x,
  condition,
  prop_na_in_condition = 0.51,
  prop_present_in_other_condition = 0.5,
  min_present_in_other_condition = 1,
  require_n_other_conditions = 1
)
```

## Arguments

- x:

  Log2 intensity matrix

- condition:

  Condition vector

- prop_na_in_condition:

  NA proportion for MNAR (default: 0.51)

- prop_present_in_other_condition:

  Required present proportion (default: 0.5)

- min_present_in_other_condition:

  Minimum present values (default: 1)

- require_n_other_conditions:

  Required conditions (default: 1)

## Value

List with keep, summary
