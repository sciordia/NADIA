# PERMANOVA R^2 via vegan::adonis2

Proportion of variance in Euclidean distances explained by the grouping.
Requires the vegan package (optional).

## Usage

``` r
.nm_permanova_r2(mat, groups, seed = 42L)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- groups:

  Factor or character vector of group labels.

- seed:

  Integer. Seed for the permutations. `adonis2()` obtains its p-value by
  permuting the group labels, so without a seed the same data give a
  different p-value on every call – and this one is exported as a column
  of
  [`nm_compute_metrics()`](https://sciordia.github.io/NADIA/reference/nm_compute_metrics.md).
  `R2` is deterministic and unaffected. The caller's RNG state is
  restored on exit.

## Value

Named list with `R2` and `p_value`, or NA if vegan unavailable.
