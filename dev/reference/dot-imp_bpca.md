# bpca: Bayesian PCA imputation (pcaMethods)

bpca: Bayesian PCA imputation (pcaMethods)

## Usage

``` r
.imp_bpca(x, args = list())
```

## Arguments

- args:

  list with optional `maxSteps` (default 100). The number of principal
  components is `ncol(x) - 1`, which `MsCoreUtils` fixes. Delegated to
  `MsCoreUtils`, which gives identical results (Bioconductor review,
  \#186); this wrapper keeps NADIA's argument names.
