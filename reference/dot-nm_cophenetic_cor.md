# MDS cophenetic correlation

Pearson correlation between the original Euclidean distances and the
distances in the 2D MDS projection. Uses scaled data (consistent with
`nm_plot_mds`).

## Usage

``` r
.nm_cophenetic_cor(mat)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

## Value

Numeric scalar (-1 to 1).
