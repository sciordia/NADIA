# NRMSE bar plot per imputation method

Bar chart of NRMSE values, ordered ascending (lower is better).

## Usage

``` r
im_plot_nrmse(metrics_df, ...)
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
im_plot_nrmse(metrics)
```
