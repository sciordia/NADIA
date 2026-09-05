# Rank normalization methods by PC1 variance explained

Computes the percentage of total variance captured by PC1 for each assay
in a SummarizedExperiment and returns a data.frame sorted in descending
order. Higher PC1 variance generally indicates stronger group
separation.

## Usage

``` r
nm_rank_pc1(
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

A `data.frame` with columns `Method`, `PC1_VarPct`, `Rank`, ordered by
`PC1_VarPct` descending.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_pc1(se, verbose = FALSE)
#>     Method PC1_VarPct Rank
#> 1 cycloess   72.57836    1
#> 2 quantile   65.48778    2
#> 3      MAD   64.86712    3
#> 4     log2   56.03859    4
```
