# Read Pattern Profiler data from a parquet file

Read Pattern Profiler data from a parquet file

## Usage

``` r
read_pattern_profiler_data(file_path, min_membership = NULL)
```

## Arguments

- file_path:

  Path to the parquet file

- min_membership:

  Optional filter by minimum membership

## Value

DataFrame with the clustering data

## Examples

``` r
if (requireNamespace("arrow", quietly = TRUE) &&
    requireNamespace("Mfuzz", quietly = TRUE) &&
    requireNamespace("e1071", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)

  # Write the clustering to parquet, then read it back
  f <- file.path(tempdir(), "pattern_profiler.parquet")
  pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                            assay_name = "Impseqrob_min",
                            auto_select_c = FALSE, c = 3,
                            output_file = f, verbose = FALSE)

  pp_data <- read_pattern_profiler_data(f, min_membership = 0.5)
  unlink(f)
  print(head(pp_data))
}
#> Filtered by membership >= 0.50: 948 -> 740 rows
#>   FeatureID Cluster Membership         A        B        D
#> 1    P14175       1   0.819835 -1.148309 0.469084 0.679225
#> 2    Q01813       1   0.819032 -1.148609 0.471726 0.676883
#> 3    P06971       1   0.818977 -1.148624 0.471853 0.676770
#> 4    P0AFI5       1   0.817308 -1.147188 0.459711 0.687478
#> 5    P63386       1   0.816561 -1.149065 0.475852 0.673212
#> 6    P64455       1   0.815585 -1.149194 0.477056 0.672138
```
