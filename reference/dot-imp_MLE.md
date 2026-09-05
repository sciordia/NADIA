# MLE: Maximum Likelihood Estimation (norm)

MLE: Maximum Likelihood Estimation (norm)

## Usage

``` r
.imp_MLE(x, args = list())
```

## Arguments

- args:

  list with optional `seed` (default 123). NAguideR strategy: no
  transpose (features x samples).

## Details

The seed is handed to
[`norm::rngseed()`](https://rdrr.io/pkg/norm/man/rngseed.html), which
seeds the internal RNG of the `norm` package rather than R's
`.Random.seed`. That state cannot be read or restored, so – unlike every
other seeded method here – this one is not wrapped in
`.rng_state()`/`.rng_restore()`. It leaves R's own RNG untouched, but
successive calls within a session are not independent.
