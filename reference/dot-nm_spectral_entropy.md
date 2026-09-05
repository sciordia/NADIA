# Spectral entropy of PCA eigenvalues

Normalized Shannon entropy of the eigenvalue distribution from PCA.
Values near 0 indicate variance concentrated in few components; values
near 1 indicate uniform spread across all components.

## Usage

``` r
.nm_spectral_entropy(mat)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

## Value

Numeric scalar (0-1).
