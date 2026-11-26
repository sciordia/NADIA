# Run differential expression analysis

Internal function that coordinates DE analysis with limma.

## Usage

``` r
.run_DE(
  se,
  comparisons,
  assay_name = NULL,
  condition_column = NULL,
  logFC = TRUE,
  logFC_up = 1,
  logFC_down = -1,
  p_adj = TRUE,
  alpha = 0.05,
  eBayes_trend = TRUE,
  eBayes_robust = TRUE,
  de_method = "limma",
  covariate_column = NULL,
  bio_replicate_column = NULL
)
```

## Arguments

- se:

  SummarizedExperiment with processed data

- comparisons:

  Comparisons to perform (result of .specify_comparisons)

- assay_name:

  Assay name to use. If NULL, uses the last available

- condition_column:

  Condition column. If NULL, uses SE metadata

- logFC:

  Apply logFC filter (default: TRUE)

- logFC_up:

  Upper logFC threshold (default: 1)

- logFC_down:

  Lower logFC threshold (default: -1)

- p_adj:

  Use adjusted p-value (default: TRUE)

- alpha:

  Significance threshold (default: 0.05)

- eBayes_trend:

  Use trend estimation in eBayes (default: TRUE)

- eBayes_robust:

  Use robust estimation in eBayes (default: TRUE)

- de_method:

  DE method: "limma" or "limpa" (default: "limma")

- covariate_column:

  Column name(s) in colData for paired/blocked design. Single string or
  character vector for multiple covariates (default: NULL)

## Value

Data frame with DE results
