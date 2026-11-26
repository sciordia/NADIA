# ACC_OI – Correlation of Original and Imputed at masked positions

Pearson correlation between the true values and the imputed values,
evaluated ONLY at the artificially masked cells (pooled across all
features and columns). Higher is better.

## Usage

``` r
.im_acc_oi(true_mat, imp_mat, na_mask)
```

## Arguments

- true_mat:

  Numeric matrix (ground truth)

- imp_mat:

  Numeric matrix (imputed)

- na_mask:

  Logical matrix

## Value

Numeric scalar (Pearson correlation over masked cells)

## Details

Note: the previous version correlated the WHOLE row (true vs imputed),
which differs in only 1-2 masked cells -\> correlation ~1 regardless of
the imputation quality, making the metric non-discriminative.
Restricting it to the masked cells measures the real accuracy.
