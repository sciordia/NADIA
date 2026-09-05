# Rank normalization methods by MDS1 variance explained

For each assay in `se`, computes the percentage of variance captured by
the first MDS dimension (using only positive eigenvalues from classical
MDS) and returns a ranking ordered by ascending MDS1 variance (lower =
better).

## Usage

``` r
nm_rank_mds1(
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

A `data.frame` with columns `Method`, `MDS1_VarPct`, `Rank`, ordered by
`MDS1_VarPct` ascending (lower = better; rank 1 = best).

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_mds1(se, verbose = FALSE)
#>     Method MDS1_VarPct Rank
#> 1 cycloess    27.18069    1
#> 2      MAD    59.70161    2
#> 3     log2    62.30919    3
#> 4 quantile    63.26223    4
```
