# Per-protein MAD averaged across groups (PRONE-style PMAD)

For each protein, computes MAD = median(\|x - median(x)\|) within each
group, then averages across groups. Returns one value per protein.

## Usage

``` r
.nm_pmad(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples)

- groups:

  Factor or character vector of group labels

## Value

Named numeric vector, one value per protein
