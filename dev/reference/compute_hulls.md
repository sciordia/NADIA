# Compute the convex hull (enclosing polygon) per group

Compute the convex hull (enclosing polygon) per group

## Usage

``` r
compute_hulls(scores_df, group_col = "Condition")
```

## Arguments

- scores_df:

  Data frame with columns PC1, PC2 and the grouping column

- group_col:

  Name of the grouping column (default: "Condition")

## Value

List of data frames, each one with columns: group, x, y
