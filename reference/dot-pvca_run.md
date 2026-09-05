# Core PVCA algorithm using lme4

1.  PCA on the expression matrix

2.  Retain PCs explaining \>= pca_threshold cumulative variance

3.  For each retained PC, fit lme4 mixed model with factors as random
    effects

4.  Extract variance components, weight by PC variance proportion

## Usage

``` r
.pvca_run(mat, annot_df, factors, pca_threshold = 0.6, verbose = TRUE)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples)

- annot_df:

  data.frame with factor columns (rows = samples)

- factors:

  Character vector of factor names

- pca_threshold:

  Numeric (0-1)

- verbose:

  Logical

## Value

data.frame with columns: label, weights
