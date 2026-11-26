# Extract results from limma fit

Extract results from limma fit

## Usage

``` r
.extract_limma_results(
  fit,
  comparisons,
  logFC_up = 1,
  logFC_down = -1,
  alpha = 0.05,
  p_adj = TRUE
)
```

## Arguments

- fit:

  limma fit object

- comparisons:

  Comparison vector

- logFC_up:

  Upper logFC threshold for "Up"

- logFC_down:

  Lower logFC threshold for "Down"

- alpha:

  Significance threshold

- p_adj:

  Use adjusted p-value (TRUE) or raw (FALSE)

## Value

Data frame with results
