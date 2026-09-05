# Hopkins statistic for clustering tendency

Manual implementation of the Hopkins statistic. Values \> 0.5 suggest
non-random clustering structure; values near 0.5 indicate uniform
distribution.

## Usage

``` r
.nm_hopkins(mat, n_sample = 10L, seed = 42L)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- n_sample:

  Integer. Number of points to sample (default 10, capped at ncol - 1).

- seed:

  Integer. Seed used for the random reference points, so the statistic
  is reproducible across calls. The caller's RNG state is restored on
  exit.

## Value

Numeric scalar (0-1).
