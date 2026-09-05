# Intra-group Pearson correlations

Computes pairwise Pearson correlations between all sample pairs within
each group, using complete observations only.

## Usage

``` r
.nm_intragroup_cor(mat, groups, method = "pearson")
```

## Arguments

- mat:

  Numeric matrix (proteins x samples)

- groups:

  Factor or character vector of group labels

## Value

Numeric vector of correlation values
