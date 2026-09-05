# softHybrid imputation: continuous sigmoid-weighted blend of MAR and MNAR

Instead of binary MAR/MNAR classification (as in combo), uses a
continuous sigmoid function to assign per-protein weights between MAR
and MNAR methods. Based on Shi et al. (bioRxiv 2026) softHybridImpute
approach.

## Usage

``` r
.impute_softHybrid(
  x,
  mar_method = "Impseqrob",
  mnar_method = "min",
  method_args = list(),
  with_value = NA_real_,
  a = 10,
  b = 5,
  lambda = 0.5,
  r0 = NULL,
  x0 = NULL
)
```

## Arguments

- x:

  Numeric matrix (proteins x samples, log2) with NAs

- mar_method:

  MAR imputation method (default: "Impseqrob", deterministic and with no
  additional dependencies; the impute_proteomics wrapper passes this
  value). Note: if a stochastic MAR/MNAR method is chosen (e.g.
  "missForest"), softHybrid does not set a seed of its own – pass one
  via method_args for reproducibility.

- mnar_method:

  MNAR imputation method (default: "min")

- method_args:

  Named list of per-method argument lists

- with_value:

  Constant value for method "with"

- a:

  Steepness of sigmoid for missing rate (default: 10)

- b:

  Steepness of sigmoid for mean intensity, in SD units (default: 5)

- lambda:

  Balance between missing rate and intensity signals (default: 0.5, in
  the range 0-1)

- r0:

  Elbow point for missing rate sigmoid (NULL = auto-detect)

- x0:

  Elbow point for intensity sigmoid (NULL = auto-detect)

## Value

List with x_imputed, weights (w_mar per protein), elbow (r0, x0),
summary
