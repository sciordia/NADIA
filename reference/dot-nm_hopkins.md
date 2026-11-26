# Hopkins statistic for clustering tendency

Manual implementation of the Hopkins statistic. Values \> 0.5 suggest
non-random clustering structure; values near 0.5 indicate uniform
distribution.

## Usage

``` r
.nm_hopkins(mat, n_sample = 10L)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- n_sample:

  Integer. Number of points to sample (default 10, capped at ncol - 1).

## Value

Numeric scalar (0-1).

## Details

The random reference points and the sampled data points are drawn from
the caller's RNG, so the statistic is reproducible only after
[`set.seed()`](https://rdrr.io/r/base/Random.html).
