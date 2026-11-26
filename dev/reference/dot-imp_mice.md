# mice: Multiple Imputation by Chained Equations

mice: Multiple Imputation by Chained Equations

## Usage

``` r
.imp_mice(x, args = list())
```

## Arguments

- args:

  list with optional `m` (default 5), `maxit` (default 5), `method`
  (default "norm") and `seed` (default `NA`). NAguideR strategy: no
  transpose, method="norm", average of m imputations.

## Details

With the default `seed = NA`,
[`mice::mice()`](https://amices.org/mice/reference/mice.html) draws from
the caller's RNG, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
before
[`impute_proteomics()`](https://sciordia.github.io/NADIA/dev/reference/impute_proteomics.md)
reproduces a run. A `seed` passed explicitly is handed to
[`mice::mice()`](https://amices.org/mice/reference/mice.html), which
then seeds the RNG itself.
