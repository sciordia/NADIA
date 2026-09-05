# Generate imputation quality metric plots and ranking

Computes NRMSE, SOR, PSS, ACC_OI via ground-truth simulation, then
generates up to 6 diagnostic plots. Each plot is wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) so that a failure
in one does not abort all.

## Usage

``` r
imputation_metrics(
  se,
  assay_name = NULL,
  methods = .IM_BENCH_METHODS,
  combo_methods = list(),
  condition_col = "Condition",
  na_prop = 0.2,
  seed = 42L,
  pattern = "random",
  method_args = list(),
  with_value = NA_real_,
  plots = "all",
  verbose = TRUE,
  output_dir = NULL,
  export_plots = TRUE,
  export_tables = TRUE,
  plot_width = 12,
  plot_height = 8,
  plot_dpi = 150
)
```

## Arguments

- se:

  SummarizedExperiment (from
  [`import_imp_matrices()`](https://sciordia.github.io/NADIA/reference/import_imp_matrices.md)
  or pipeline).

- assay_name:

  Character. Assay name to use as starting point. NULL (default) = first
  assay.

- methods:

  Character vector of individual methods to benchmark. Default:
  `.IM_BENCH_METHODS` (16 methods). Use `NULL` or `character(0)` to skip
  individual methods when only using combo_methods.

- combo_methods:

  Named list of combo/softHybrid configurations. Each element:
  `list(mar_method, mnar_method, mode)`. See
  [`im_compute_metrics()`](https://sciordia.github.io/NADIA/reference/im_compute_metrics.md)
  for details.

- condition_col:

  Character. Column in colData for combo methods. Default `"Condition"`.

- na_prop:

  Numeric (0-1). Proportion of artificial NAs. Default 0.20. Only used
  when `pattern = "random"`. Ignored when `pattern = "from_data"` (NA
  proportion is derived from the data, replicating NAguideR).

- seed:

  Integer. Random seed. Default 42.

- pattern:

  Character. `"random"` or `"from_data"`. Default `"random"`.

- method_args:

  Named list of per-method argument lists.

- with_value:

  Constant for method `"with"`.

- plots:

  Character vector of plot names or `"all"` (default). Valid: `"nrmse"`,
  `"sor"`, `"pss"`, `"acc_oi"`, `"ranking"`, `"metrics"`.

- verbose:

  Logical. Print progress. Default TRUE.

- output_dir:

  Character. Path to export directory. If `NULL` (default), no files are
  exported. When set, tables (TSV) and/or plots (PNG) are saved to this
  directory (created if needed).

- export_plots:

  Logical. Export plots as PNG when `output_dir` is set. Default `TRUE`.

- export_tables:

  Logical. Export tables as TSV when `output_dir` is set. Default
  `TRUE`.

- plot_width:

  Numeric. Width in inches for exported plots. Default `12`.

- plot_height:

  Numeric. Height in inches for exported plots. Default `8`.

- plot_dpi:

  Numeric. Resolution for exported plots. Default `150`.

## Value

Named list of ggplot objects (or NULL for failed plots), plus
`metrics_table`: a data.frame from
[`im_compute_metrics()`](https://sciordia.github.io/NADIA/reference/im_compute_metrics.md).

## Examples

``` r
data(nadia_dia)
se <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)

# The default sweeps 16 methods; a handful is enough to illustrate it
res <- imputation_metrics(se, methods = c("min", "zero", "knn"),
                          plots = "ranking", verbose = FALSE)
res$metrics_table
#>   Method     NRMSE  SOR         PSS    ACC_OI NRMSE_Rank SOR_Rank PSS_Rank
#> 1    knn 0.1293948 1349 0.004570796 0.9916164          1        1        1
#> 2    min 3.1934374 2698 0.456162758        NA          2        2        2
#> 3   zero 6.3558709 4047 0.470781473        NA          3        3        3
#>   ACC_OI_Rank Rank_Mean
#> 1           1         1
#> 2          NA         2
#> 3          NA         3

# Two-stage combos can be compared alongside the individual methods
res2 <- imputation_metrics(se, methods = c("min", "knn"),
  combo_methods = list(
    "Impseq+min" = list(mar_method = "Impseq", mnar_method = "min")),
  plots = "ranking", verbose = FALSE)
res2$metrics_table$Method
#> [1] "knn"        "Impseq+min" "min"       
```
