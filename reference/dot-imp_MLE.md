# MLE: Maximum Likelihood Estimation (norm)

MLE: Maximum Likelihood Estimation (norm)

## Usage

``` r
.imp_MLE(x, args = list())
```

## Arguments

- args:

  list with optional `seed` (default `NULL`). NAguideR strategy: no
  transpose (features x samples).

## Details

`norm` has an RNG of its own, separate from R's, which
[`norm::rngseed()`](https://rdrr.io/pkg/norm/man/rngseed.html)
initialises. When no `seed` is given, it is drawn from the caller's RNG,
so [`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`impute_proteomics()`](https://sciordia.github.io/NADIA/reference/impute_proteomics.md)
reproduces a run here as it does for every other stochastic method. A
`seed` passed explicitly is handed to
[`norm::rngseed()`](https://rdrr.io/pkg/norm/man/rngseed.html)
unchanged. One limitation belongs to `norm` and no seed removes it: its
generator keeps state that `rngseed()` does not reset, so a second `MLE`
imputation in the same R session differs from the first even with the
same seed. A script run in a fresh session is reproducible.
