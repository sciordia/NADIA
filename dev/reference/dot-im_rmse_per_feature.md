# Per-feature RMSE for SOR computation

For each row (feature) with artificial NAs, computes RMSE between true
and imputed values. Returns a named numeric vector.

## Usage

``` r
.im_rmse_per_feature(true_mat, imp_mat, na_mask)
```

## Arguments

- true_mat:

  Numeric matrix (ground truth)

- imp_mat:

  Numeric matrix (imputed)

- na_mask:

  Logical matrix

## Value

Named numeric vector of per-feature RMSE
