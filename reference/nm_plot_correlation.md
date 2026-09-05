# Intra-group correlation boxplot per method (PRONE-style)

Pairwise within-group correlations for each method shown as boxplots
with error bars. Higher and less variable values indicate better
normalization.

## Usage

``` r
nm_plot_correlation(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  cor_method = "pearson",
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

- cor_method:

  Correlation method passed to
  [`cor()`](https://rdrr.io/r/stats/cor.html): `"pearson"`,
  `"spearman"`, or `"kendall"`. Default `"pearson"`.

- ...:

  Additional arguments (unused).

## Value

ggplot object.

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                            verbose = FALSE)

nm_plot_correlation(se, cor_method = "pearson")

```
