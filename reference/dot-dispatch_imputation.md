# Central dispatcher for imputation methods

Central dispatcher for imputation methods

## Usage

``` r
.dispatch_imputation(x, method, method_args = list(), with_value = NA_real_)
```

## Arguments

- x:

  Numeric matrix (proteins x samples)

- method:

  Method name (one of .IMP_METHODS_ALL minus "combo")

- method_args:

  Named list of per-method argument lists

- with_value:

  Constant value for method "with"

## Value

Imputed matrix
