# Evaluate a range of cluster numbers

Evaluate a range of cluster numbers

## Usage

``` r
evaluate_cluster_range(
  eset_std,
  c_range = 2:10,
  m = NULL,
  n_starts = 3L,
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

- n_starts:

  Number of random starts per candidate number of clusters; the metrics
  are averaged over them. Drawn from the caller's RNG.

- verbose:

  Show progress

## Value

DataFrame with the metrics per number of clusters
