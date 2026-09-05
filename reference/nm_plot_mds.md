# MDS 2D scatter per method

Classical MDS on Euclidean distances of scaled, complete-case data.
Faceted by method, colored by condition.

## Usage

``` r
nm_plot_mds(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  mds_scales = c("free", "fixed"),
  ...
)
```

## Arguments

- se:

  SummarizedExperiment from
  [`import_norm_matrices()`](https://sciordia.github.io/NADIA/reference/import_norm_matrices.md).

- assay_names:

  Character vector of assay names to include. NULL = all.

- condition_col:

  Column in colData with condition labels.

- mds_scales:

  Facet scaling: `"free"` (default) allows independent axes per panel;
  `"fixed"` forces shared axes for easier cross-method comparison.

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

nm_plot_mds(se)

```
