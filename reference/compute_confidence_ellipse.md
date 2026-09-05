# Compute the confidence ellipse per group

Computes the coordinates of a confidence ellipse based on the
chi-squared distribution, similar to FactoMineR::coord.ellipse and
ggplot2::stat_ellipse.

## Usage

``` r
compute_confidence_ellipse(
  scores_df,
  group_col = "Condition",
  level = 0.95,
  npoints = 100
)
```

## Arguments

- scores_df:

  Data frame with columns PC1, PC2 and the grouping column

- group_col:

  Name of the grouping column (default: "Condition")

- level:

  Confidence level (default: 0.95)

- npoints:

  Number of points used to draw the ellipse (default: 100)

## Value

List of data frames, each one with columns: group, x, y
