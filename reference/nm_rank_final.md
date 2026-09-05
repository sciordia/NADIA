# Compute a combined final ranking across PCV, PMAD, PEV, Correlation and group-separation (PC1 F-ratio)

For each method, the final rank is the mean of five individual ranks:
four intragroup-precision metrics (pcv, pmad, pev, cor) plus a
group-separation metric (`Rank_Sep`, the ANOVA F-ratio on PC1 scores,
higher = better separation). Lower final rank = better overall
normalization quality.

## Usage

``` r
nm_rank_final(
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

  Correlation method forwarded to
  [`nm_rank_cor()`](https://sciordia.github.io/NADIA/reference/nm_rank_cor.md).

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `data.frame` with columns: `Method`, `Rank_PCV`, `Rank_PMAD`,
`Rank_PEV`, `Rank_Cor`, `Rank_Sep`, `Rank_Final`, ordered by
`Rank_Final` ascending (best first).

## Details

Note: earlier versions averaged `Rank_PC1` (PC1 variance %, higher =
better) and `Rank_MDS1` (MDS1 variance %, lower = better). Both capture
the same dominant axis of variation but with opposite directions, so
they partially cancelled within the mean. They are replaced here by a
single, unambiguous group-separation rank based on the PC1 F-ratio. The
standalone
[`nm_rank_pc1()`](https://sciordia.github.io/NADIA/reference/nm_rank_pc1.md)
/
[`nm_rank_mds1()`](https://sciordia.github.io/NADIA/reference/nm_rank_mds1.md)
functions are unchanged for individual inspection.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_rank_final(se, verbose = FALSE)
#>     Method Rank_PCV Rank_PMAD Rank_PEV Rank_Cor Rank_Sep Rank_Final
#> 1 cycloess        1         1        1        1        3        1.4
#> 2 quantile        3         2        2        4        1        2.4
#> 3      MAD        2         3        3        3        2        2.6
#> 4     log2        4         4        4        2        4        3.6
```
