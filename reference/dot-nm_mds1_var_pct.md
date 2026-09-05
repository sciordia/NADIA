# Percentage of variance explained by MDS dimension 1

Computes the percentage of variance captured by the first MDS dimension,
using only positive eigenvalues from classical MDS
(`cmdscale(eig = TRUE)`). Uses scaled data (consistent with
`nm_plot_mds`).

## Usage

``` r
.nm_mds1_var_pct(mat)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

## Value

Numeric scalar (0-100).
