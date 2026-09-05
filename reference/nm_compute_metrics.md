# Compute quantitative group-separation metrics per normalization method

Iterates over assays in a SummarizedExperiment and computes ten metrics
that quantify how well the normalization separates sample groups.

## Usage

``` r
nm_compute_metrics(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  seed = 42L,
  ...
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

- seed:

  Integer. Seed for the two metrics that use randomness: the Hopkins
  statistic, which samples random reference points, and the PERMANOVA
  p-value, which permutes the group labels. Without it neither column is
  reproducible between runs. The caller's RNG state is restored on exit.
  Default 42.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row per method and columns: `Method`,
`PC1_VarPct`, `PC1_F_ratio`, `PERMANOVA_R2`, `PERMANOVA_pval`,
`Silhouette_mean`, `MDS_GOF`, `MDS_CophCor`, `Spectral_Entropy`,
`CumVar_PC2`, `Hopkins`, `Condition_Number`.

## Details

- **PC1_VarPct**: % variance explained by PC1 (higher = more structure).

- **PC1_F_ratio**: ANOVA F-ratio on PC1 scores (higher = better
  separation).

- **PERMANOVA_R2**: Proportion of variance explained by grouping
  (requires vegan; NA if not installed).

- **Silhouette_mean**: Mean silhouette width (requires cluster; NA if
  not installed).

- **MDS_GOF**: Goodness-of-fit of 2D MDS projection. Note: in
  benchmarking context, *lower* values tend to indicate better
  normalization (a good method preserves multi-dimensional biological
  variability that 2D cannot capture).

- **MDS_CophCor**: Cophenetic correlation between original and MDS
  distances. Same caveat as MDS_GOF: *lower* values tend to correlate
  with better benchmarking performance.

- **Spectral_Entropy**: Normalized Shannon entropy of PCA eigenvalues
  (0-1). Values near 0 = variance concentrated in few PCs; near 1 =
  uniform spread.

- **CumVar_PC2**: Cumulative % variance explained by PC1 + PC2 (0-100).

- **Hopkins**: Hopkins statistic for clustering tendency (0-1). Values
  \> 0.5 suggest non-random cluster structure.

- **Condition_Number**: Ratio of largest to smallest PCA eigenvalue (\>=
  1). High values indicate multicollinearity or numerical instability.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

metrics_df <- nm_compute_metrics(se)
metrics_df[, c("Method", "PC1_VarPct", "PC1_F_ratio", "Silhouette_mean")]
#>     Method PC1_VarPct PC1_F_ratio Silhouette_mean
#> 1     log2   56.03859    120.9196       0.2118819
#> 2 cycloess   72.57836   9790.2845       0.3790455
#> 3 quantile   65.48778  10168.4649       0.3760405
#> 4      MAD   64.86712   9803.3024       0.3862026
```
