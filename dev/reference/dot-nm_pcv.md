# Per-protein CV averaged across groups (PRONE-style PCV)

For each protein, computes CV = 100 \* SD / \|mean\| within each group,
then averages across groups. Returns one value per protein.

## Usage

``` r
.nm_pcv(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), log2 scale

- groups:

  Factor or character vector of group labels (length = ncol(mat))

## Value

Named numeric vector, one value per protein
