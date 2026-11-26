# PERMANOVA R^2 via vegan::adonis2

Proportion of variance in Euclidean distances explained by the grouping.
Requires the vegan package (optional).

## Usage

``` r
.nm_permanova_r2(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- groups:

  Factor or character vector of group labels.

## Value

Named list with `R2` and `p_value`, or NA if vegan unavailable.

## Details

`adonis2()` obtains its p-value by permuting the group labels with the
caller's RNG, so the p-value is reproducible only after
[`set.seed()`](https://rdrr.io/r/base/Random.html). `R2` is
deterministic.
