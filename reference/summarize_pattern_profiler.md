# Summary of the Pattern Profiler data

Summary of the Pattern Profiler data

## Usage

``` r
summarize_pattern_profiler(data)
```

## Arguments

- data:

  Pattern Profiler DataFrame

## Value

List with summary statistics

## Examples

``` r
if (requireNamespace("Mfuzz", quietly = TRUE) &&
    requireNamespace("e1071", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)
  pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                  assay_name = "Impseqrob_min",
                                  auto_select_c = FALSE, c = 3,
                                  verbose = FALSE)

  info <- summarize_pattern_profiler(pp$long_output)
  print(info$cluster_summary)
}
#> # A tibble: 3 × 6
#>   Cluster n_entries n_unique_features mean_membership min_membership
#>     <int>     <int>             <int>           <dbl>          <dbl>
#> 1       1       412               412           0.860          0.250
#> 2       2       177               177           0.476          0.252
#> 3       3       359               359           0.707          0.255
#> # ℹ 1 more variable: max_membership <dbl>
```
