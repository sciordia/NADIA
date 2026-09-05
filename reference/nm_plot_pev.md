# PEV boxplot per method (PRONE-style)

One boxplot per normalization method showing the per-protein variance
distribution (variance averaged across groups for each protein). Lower
values indicate better within-group consistency. With `diff = TRUE`,
shows % reduction vs `baseline` as a bar chart.

## Usage

``` r
nm_plot_pev(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  diff = FALSE,
  baseline = "log2",
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

- diff:

  Logical. If TRUE, show % reduction vs `baseline`. Default FALSE.

- baseline:

  Character. Assay name used as reference for diff mode. Default
  `"log2"`.

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

nm_plot_pev(se)

```
