# Build the output table in LONG format

Creates one row per FeatureID-Cluster combination where membership \>=
min_membership. This lets a protein appear in several clusters (soft
clustering).

## Usage

``` r
build_long_output(cl, eset_std, conditions, min_membership)
```

## Arguments

- cl:

  Mfuzz clustering object

- eset_std:

  Standardised ExpressionSet (with z-scores)

- conditions:

  Vector of condition names

- min_membership:

  Minimum membership threshold

## Value

DataFrame in long format
