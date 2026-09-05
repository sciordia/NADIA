# Re-impute a matrix with multiple methods

Iterates over a vector of imputation method names, calling
[`.dispatch_imputation()`](https://sciordia.github.io/NADIA/reference/dot-dispatch_imputation.md)
from Imputation.R for each. Also supports combo and softHybrid methods
via `combo_methods` parameter. Methods that fail are omitted with a
warning.

## Usage

``` r
.im_reimpute(
  mat_with_na,
  methods = character(0),
  combo_methods = list(),
  condition = NULL,
  method_args = list(),
  with_value = NA_real_,
  verbose = TRUE
)
```

## Arguments

- mat_with_na:

  Numeric matrix with artificial NAs

- methods:

  Character vector of individual method names

- combo_methods:

  Named list of combo/softHybrid configurations. Each element is a list
  with: `mode` ("combo" or "softHybrid", default "combo"), `mar_method`,
  `mnar_method`, and optionally other parameters. The list name is used
  as the method label.

- condition:

  Factor or character vector of conditions (required for combo methods,
  ignored for individual methods).

- method_args:

  Named list of per-method argument lists

- with_value:

  Constant for method "with"

- verbose:

  Logical. Print progress. Default TRUE.

## Value

Named list of imputed matrices (one per successful method)
