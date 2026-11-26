# Build an extended combined data.frame from opdea + benchmark metrics

Merges OpDEA metrics and benchmark classification metrics by Assay +
Comparison, converts Performance to numeric, and selects the 11 extended
metrics for ranking.

## Usage

``` r
.bm_build_extended(opdea_combined, bench_metrics_combined)
```

## Arguments

- opdea_combined:

  data.frame from import_opdea_results()

- bench_metrics_combined:

  data.frame from import_benchmark_metrics()

## Value

data.frame with columns: Assay, Comparison, and the 11 metrics in
.BM_EXTENDED_METRICS
