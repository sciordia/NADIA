# Average silhouette width

Mean silhouette width when samples are clustered by condition labels.
Requires the cluster package (optional).

## Usage

``` r
.nm_silhouette_avg(mat, groups)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs.

- groups:

  Factor or character vector of group labels.

## Value

Numeric scalar (-1 to 1), or NA if cluster unavailable.
