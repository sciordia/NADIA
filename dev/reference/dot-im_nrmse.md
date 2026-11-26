# NRMSE – Normalized Root Mean Squared Error

Computes NRMSE between true and imputed values at positions indicated by
the NA mask. Lower is better.

## Usage

``` r
.im_nrmse(true_mat, imp_mat, na_mask)
```

## Arguments

- true_mat:

  Numeric matrix (ground truth, complete)

- imp_mat:

  Numeric matrix (imputed from artificially-NA'd version)

- na_mask:

  Logical matrix (TRUE = position was set to NA)

## Value

Numeric scalar (NRMSE)
