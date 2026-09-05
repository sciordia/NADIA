# PCA scatter plot per method

PC1 vs PC2 scatter, colored by condition, faceted by method. Percentage
of variance explained shown on each axis.

## Usage

``` r
nm_plot_pca(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  pca_scales = c("free", "fixed"),
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

- pca_scales:

  Facet scaling: `"free"` (default) allows independent axes per method;
  `"fixed"` uses shared axes to compare separation magnitude.

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

nm_plot_pca(se, pca_scales = "fixed")

```
