# Horizontal bar chart of PC1 variance ranking

Produces a horizontal bar chart with normalization methods ordered by
descending PC1 variance percentage. Uses the same PRONE-style palette as
other `nm_plot_*()` functions.

## Usage

``` r
nm_plot_pc1_ranking(se, assay_names = NULL, condition_col = "Condition", ...)
```

## Arguments

- se:

  SummarizedExperiment from
  [`import_norm_matrices()`](https://sciordia.github.io/NADIA/reference/import_norm_matrices.md).

- assay_names:

  Character vector of assay names to include. NULL = all.

- condition_col:

  Column in colData with condition labels.

- ...:

  Additional arguments (unused).

## Value

ggplot object.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_plot_pc1_ranking(se)

```
