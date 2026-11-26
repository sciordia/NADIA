# Prepare benchmark data

Filters de_res by assay and comparisons, validates inputs, resolves
p-value column.

## Usage

``` r
.prepare_benchmark_data(
  de_res,
  expected_values,
  alpha = 0.05,
  lfc_thr = 0,
  p_col = "adj.P.Val",
  comparisons = NULL,
  assay = NULL,
  species_df = NULL
)
```

## Arguments

- de_res:

  Data frame with DE results (must include Species column)

- expected_values:

  Data frame with expected values

- alpha:

  Significance threshold

- lfc_thr:

  Log fold-change threshold

- p_col:

  P-value column name

- comparisons:

  Comparisons to include (NULL = all)

- assay:

  Assay to filter (NULL = all)

- species_df:

  Data frame with Protein.IDs and Species columns (optional). If
  provided, Species is merged into de_res. If NULL, de_res must already
  contain a Species column.

## Value

List with filtered de_res, expected_values, and resolved p_col
