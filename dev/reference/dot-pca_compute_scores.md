# Compute PCA scores from a SE assay

Compute PCA scores from a SE assay

## Usage

``` r
.pca_compute_scores(
  se,
  assay_name,
  na_action = "complete",
  fill_value = -1,
  center = TRUE,
  scale. = TRUE,
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment

- assay_name:

  Character. Assay to use.

- na_action:

  Character: "complete" (remove rows with NAs) or "fill".

- fill_value:

  Numeric. Replacement for NAs when na_action="fill".

- center:

  Logical. Center variables before PCA (default TRUE).

- scale.:

  Logical. Scale variables before PCA (default TRUE).

- verbose:

  Logical.

## Value

List with: scores (data.frame PC1, PC2, Sample + colData), pct_var
(numeric vector of % variance per PC), n_proteins (integer).
