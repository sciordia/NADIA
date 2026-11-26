# Per-protein variance averaged across groups (PRONE-style PEV)

For each protein, computes variance within each group, then averages
across groups. Returns one value per protein.

## Usage

``` r
.nm_pev(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples)

- groups:

  Factor or character vector of group labels

## Value

Named numeric vector, one value per protein
