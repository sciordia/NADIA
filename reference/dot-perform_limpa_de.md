# Run limpa dpcDE analysis

Uses limpa::dpcDE with precision weights from the EList (standard
errors) produced by dpcQuantByRow. The EList is stored in SE metadata by
impute_proteomics() when imp_method="limpa".

## Usage

``` r
.perform_limpa_de(
  elist,
  condition_vector,
  comparisons,
  covariate = NULL,
  block = NULL,
  eBayes_trend = FALSE,
  eBayes_robust = FALSE
)
```

## Arguments

- elist:

  EList from limpa (with \$E and \$weights)

- condition_vector:

  Condition vector aligned with columns

- comparisons:

  Comparison vector

- covariate:

  Optional covariate for the model

- eBayes_trend:

  Use trend estimation in eBayes (default: FALSE, vooma already models
  the trend)

- eBayes_robust:

  Use robust estimation in eBayes (default: FALSE,
  voomaLmFitWithImputation already handles imputed proteins)

## Value

limma MArrayLM fit object
