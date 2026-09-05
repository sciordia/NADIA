# Rank normalization methods by median PMAD (ascending – lower is better)

Rank normalization methods by median PMAD (ascending – lower is better)

## Usage

``` r
nm_rank_pmad(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  verbose = TRUE
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

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `data.frame` with columns `Method`, `Median_PMAD`, `Rank`.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_pmad(se, verbose = FALSE)
#>     Method Median_PMAD Rank
#> 1 cycloess  0.06634680    1
#> 2 quantile  0.06851739    2
#> 3      MAD  0.07282100    3
#> 4     log2  0.17211589    4
```
