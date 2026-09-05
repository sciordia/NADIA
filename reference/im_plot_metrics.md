# Faceted bar chart of all 4 imputation quality metrics

Similar to
[`nm_plot_metrics()`](https://sciordia.github.io/NADIA/reference/nm_plot_metrics.md):
one facet per metric with free y-scales, labeled values on bars.

## Usage

``` r
im_plot_metrics(metrics_df, ...)
```

## Arguments

- metrics_df:

  data.frame from
  [`im_compute_metrics()`](https://sciordia.github.io/NADIA/reference/im_compute_metrics.md).

- ...:

  Additional arguments (unused).

## Value

ggplot object.

## Examples

``` r
data(nadia_dia)
se <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)
metrics <- im_compute_metrics(se, methods = c("min", "zero", "knn"),
                              verbose = FALSE)
im_plot_metrics(metrics)
```
