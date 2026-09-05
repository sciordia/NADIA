# PSS – Procrustes Statistical Shape analysis

PCA on both true and imputed matrices (transposed: samples as rows),
retain components up to 95% cumulative variance, then Procrustes SS.
Requires vegan (optional). Lower is better.

## Usage

``` r
.im_pss(true_mat, imp_mat)
```

## Arguments

- true_mat:

  Numeric matrix (ground truth, proteins x samples)

- imp_mat:

  Numeric matrix (imputed)

## Value

Numeric scalar (Procrustes SS), or NA if vegan not available
