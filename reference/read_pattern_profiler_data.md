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
#>   FeatureID Cluster Membership        A        B         D
#> 1    P32610       1   0.977237 0.638391 0.514077 -1.152468
#> 2    Q04182       1   0.977214 0.638649 0.513799 -1.152448
#> 3    P21375       1   0.977183 0.638171 0.514313 -1.152484
#> 4    P36090       1   0.977170 0.638774 0.513665 -1.152439
#> 5    P21147       1   0.977150 0.638096 0.514394 -1.152490
#> 6    Q01476       1   0.977112 0.638027 0.514468 -1.152495
```
