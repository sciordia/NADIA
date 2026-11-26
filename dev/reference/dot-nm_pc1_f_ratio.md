# PC1 F-ratio (between-group / within-group variance on PC1 scores)

One-way ANOVA F-statistic on PC1 scores grouped by condition. Higher
values indicate better group separation along PC1.

## Usage

``` r
.nm_pc1_f_ratio(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- groups:

  Factor or character vector of group labels.

## Value

Numeric scalar (0-Inf).
