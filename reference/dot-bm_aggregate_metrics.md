# Aggregate metrics by Assay

Aggregate metrics by Assay

## Usage

``` r
.bm_aggregate_metrics(df, metrics, agg_fun = mean)
```

## Arguments

- df:

  Combined opdea data.frame (with Assay column)

- metrics:

  Character vector of metric column names

- agg_fun:

  Aggregation function (mean or median)

## Value

data.frame with one row per Assay
