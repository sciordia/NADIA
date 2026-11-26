# Condition number of the covariance matrix (via PCA eigenvalues)

Ratio of the largest to smallest non-zero eigenvalue from PCA of the
sample covariance matrix. High values indicate numerical instability or
multicollinearity.

## Usage

``` r
.nm_condition_number(mat)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

## Value

Numeric scalar \>= 1 (Inf if smallest eigenvalue is zero).
