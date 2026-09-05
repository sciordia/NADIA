# Compute imputation quality metrics (NRMSE, SOR, PSS, ACC_OI)

Extracts complete-case rows from the assay, introduces artificial NAs,
re-imputes with each method, and computes four quality metrics plus a
combined ranking.

## Usage

``` r
im_compute_metrics(
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
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment (from
  [`import_imp_matrices()`](https://sciordia.github.io/NADIA/reference/import_imp_matrices.md)
  or pipeline)

- assay_name:

  Character. Name of the assay to use as starting point (typically the
  normalized assay, pre-imputation, with real NAs).

- methods:

  Character vector of individual imputation methods to benchmark.
  Default: `.IM_BENCH_METHODS` (16 methods). Use `character(0)` or
  `NULL` to skip individual methods.

- combo_methods:

  Named list of combo/softHybrid configurations. Each element is a list
  with: `mode` ("combo" or "softHybrid", default "combo"), `mar_method`,
  `mnar_method`. The list name is used as the method label in results.
  Default: empty list (no combo methods).

- condition_col:

  Character. Column in colData with condition labels, required for combo
  methods. Default `"Condition"`.

- na_prop:

  Numeric (0-1). Proportion of artificial NAs. Default 0.20. Only used
  when `pattern = "random"`. Ignored when `pattern = "from_data"` (NA
  proportion is derived from the data, replicating NAguideR).

- seed:

  Integer. Random seed. Default 42.

- pattern:

  Character. NA introduction pattern: `"random"` for uniform random NAs,
  or `"from_data"` to replicate NAguideR's simulation strategy (row
  ratio + per-column distribution from the actual data).

- method_args:

  Named list of per-method argument lists.

- with_value:

  Constant value for method `"with"`.

- verbose:

  Logical. Print progress. Default TRUE.

## Value

data.frame with columns: Method, NRMSE, SOR, PSS, ACC_OI, NRMSE_Rank,
SOR_Rank, PSS_Rank, ACC_OI_Rank, Rank_Mean. Ordered by Rank_Mean (best
first).

## Examples

``` r
data(nadia_dia)
se <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)

# Ground-truth simulation: complete rows are masked and re-imputed, so the
# error of each method can be measured against the values it did not see
metrics <- im_compute_metrics(se, methods = c("min", "zero", "knn"),
                              verbose = FALSE)
metrics[, c("Method", "NRMSE", "SOR", "Rank_Mean")]
#>   Method     NRMSE  SOR Rank_Mean
#> 1    knn 0.1293948 1349         1
#> 2    min 3.1934374 2698         2
#> 3   zero 6.3558709 4047         3
```
