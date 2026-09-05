# missForest: Random Forest imputation

missForest: Random Forest imputation

## Usage

``` r
.imp_missForest(x, args = list())
```

## Arguments

- args:

  list with optional `maxiter` (default 10), `ntree` (default 100),
  `mtry` (default floor(nrow(x)^(1/3)), NAguideR cube root strategy)

## Details

This method is stochastic and takes no seed of its own, because
[`missForest::missForest()`](https://rdrr.io/pkg/missForest/man/missForest.html)
does not accept one. To reproduce a run, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`impute_proteomics()`](https://sciordia.github.io/NADIA/reference/impute_proteomics.md);
unlike the seeded methods, this one draws from – and therefore advances
– the caller's RNG.
