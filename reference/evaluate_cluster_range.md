# Evaluate a range of cluster numbers

Evaluate a range of cluster numbers

## Usage

``` r
evaluate_cluster_range(
  eset_std,
  c_range = 2:10,
  m = NULL,
  seeds = c(42, 123, 456),
  verbose = TRUE
)
```

## Arguments

- eset_std:

  Standardised ExpressionSet

- c_range:

  Vector of cluster numbers to evaluate

- m:

  Fuzziness parameter

- seeds:

  Seeds for reproducibility

- verbose:

  Show progress

## Value

DataFrame with the metrics per number of clusters
