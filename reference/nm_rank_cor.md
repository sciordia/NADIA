# Rank normalization methods by median intragroup correlation (descending – higher is better)

Rank normalization methods by median intragroup correlation (descending
– higher is better)

## Usage

``` r
nm_rank_cor(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  cor_method = "pearson",
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

- cor_method:

  Correlation method: "pearson", "spearman", or "kendall".

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `data.frame` with columns `Method`, `Median_Cor`, `Rank`.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_cor(se, cor_method = "pearson", verbose = FALSE)
#>     Method Median_Cor Rank
#> 1 cycloess  0.9866091    1
#> 2     log2  0.9865677    2
#> 3      MAD  0.9865677    3
#> 4 quantile  0.9864951    4
```
