# Check which interaction terms are safe (n_levels \< n_obs)

Returns only interactions whose combined factor levels are strictly less
than the number of observations. This prevents lme4 from failing on
saturated random effects (e.g., Condition:Patient in paired designs).

## Usage

``` r
.pvca_safe_interactions(annot_df, factors, verbose = TRUE)
```

## Arguments

- annot_df:

  data.frame with factor columns

- factors:

  Character vector of main effect names

- verbose:

  Logical

## Value

Character vector of safe interaction terms (e.g., "A:B")
