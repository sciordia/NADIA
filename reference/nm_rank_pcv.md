# Rank normalization methods by median PCV (ascending – lower is better)

Rank normalization methods by median PCV (ascending – lower is better)

## Usage

``` r
nm_rank_pcv(
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

A `data.frame` with columns `Method`, `Median_PCV`, `Rank`.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_pcv(se, verbose = FALSE)
#>     Method Median_PCV Rank
#> 1 cycloess  0.9442522    1
#> 2      MAD  0.9623733    2
#> 3 quantile  0.9686333    3
#> 4     log2  1.9240485    4
```
