# Run limma analysis

Run limma analysis

## Usage

``` r
.perform_limma(
  data,
  condition_vector,
  comparisons,
  covariate = NULL,
  block = NULL,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE
)
```

## Arguments

- data:

  Log2 intensity matrix

- condition_vector:

  Condition vector aligned with columns

- comparisons:

  Comparison vector

- covariate:

  Optional covariate for the model

- eBayes_trend:

  Use trend estimation in eBayes (default: TRUE)

- eBayes_robust:

  Use robust estimation in eBayes (default: TRUE)

## Value

limma fit object
