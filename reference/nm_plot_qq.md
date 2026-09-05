# Q-Q plot per method (NormalyzerDE style)

Quantile-quantile plot against a normal distribution for a single
sample, one facet per normalization method. Default uses the first
column (as NormalyzerDE does).

## Usage

``` r
nm_plot_qq(
  se,
  assay_names = NULL,
  condition_col = "Condition",
  which_sample = NULL,
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

- which_sample:

  Character. Name of the sample to inspect. Default = first column of
  `se`.

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

nm_plot_qq(se, which_sample = "A_1")

```
