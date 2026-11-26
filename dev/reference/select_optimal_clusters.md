# Select the optimal number of clusters

Select the optimal number of clusters

## Usage

``` r
select_optimal_clusters(
  eset_std,
  c_range = 2:10,
  m = NULL,
  method = c("xb", "consensus", "elbow"),
  n_starts = 3L,
  verbose = TRUE
)
```

## Arguments

- eset_std:

  Standardised ExpressionSet

- c_range:

  Range of cluster numbers to evaluate

- m:

  Fuzziness parameter

- method:

  Method: "xb", "consensus", "elbow"

- n_starts:

  Random starts per candidate, passed on to
  [`evaluate_cluster_range()`](https://sciordia.github.io/NADIA/dev/reference/evaluate_cluster_range.md)

- verbose:

  Show progress

## Value

List with optimal_c, metrics, m
