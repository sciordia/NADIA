# PI: Perseus-style imputation (down-shifted normal distribution)

PI: Perseus-style imputation (down-shifted normal distribution)

## Usage

``` r
.imp_PI(x, args = list())
```

## Arguments

- args:

  list with optional `width` (default 0.3) and `downshift` (default
  1.8).

## Details

The imputed values are drawn with
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) from the caller's RNG.
To reproduce a run, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`impute_proteomics()`](https://sciordia.github.io/NADIA/reference/impute_proteomics.md).
