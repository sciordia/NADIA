# MDS goodness-of-fit

GOF\[1\] from [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) with
`eig = TRUE`: proportion of variance retained in the 2D MDS projection.
Uses scaled data (consistent with `nm_plot_mds`).

## Usage

``` r
.nm_mds_gof(mat)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

## Value

Numeric scalar (0-1).
