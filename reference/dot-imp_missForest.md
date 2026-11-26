# missForest: Random Forest imputation

Which unit the forests learn from is set by `observations`:

- `"features"` (default): the proteins are the observations and the
  samples the variables, so each sample is predicted from the other
  samples, with forests trained over thousands of proteins. This is the
  default of
  [`MsCoreUtils::impute_RF()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  and the orientation of `knn` and `bpca` here.

- `"samples"`: the samples are the observations and every protein is a
  variable, as NAguideR does. Each forest then trains on at most n - 1
  samples, which suits large cohorts only.

## Usage

``` r
.imp_missForest(x, args = list())
```

## Arguments

- args:

  list with optional `observations` (`"features"`, the default, or
  `"samples"`), `maxiter` (default 10), `ntree` (default 100) and `mtry`
  (default `floor(sqrt(p))` with `"features"`, the missForest default,
  and `floor(p^(1/3))` with `"samples"`, NAguideR's cube root of the
  number of proteins). The fit is delegated to
  [`MsCoreUtils::impute_RF()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html),
  which gives identical results.

## Details

Until NADIA 0.99.2 the only orientation was `"samples"`. On `nadia_dia`
(12 samples) with values hidden at random, the RMSE was 0.264 against
0.402, and on the TMT example 0.113 against 0.165 with 24 samples and
0.153 against 0.455 with 12 (Bioconductor review, \#186).

This method is stochastic and takes no seed of its own, because
[`missForest::missForest()`](https://rdrr.io/pkg/missForest/man/missForest.html)
does not accept one. To reproduce a run, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`impute_proteomics()`](https://sciordia.github.io/NADIA/reference/impute_proteomics.md),
as for every stochastic method in NADIA.
