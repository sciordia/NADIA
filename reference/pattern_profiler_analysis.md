# Pattern Profiler Analysis

Complete clustering pipeline starting from a SummarizedExperiment.
Writes a parquet file in LONG format for visualisation.

## Usage

``` r
pattern_profiler_analysis(
  se_proc,
  DEPs_results,
  assay_name = NULL,
  filter_mode = c("any", "all", "specific"),
  alpha = 0.05,
  comparison = NULL,
  condition_order = NULL,
  aggregate = c("median", "mean"),
  c_range = 2:10,
  auto_select_c = TRUE,
  c = NULL,
  selection_method = c("xb", "consensus", "elbow"),
  min_membership = 0.25,
  seed = 42,
  seeds = c(42, 123, 456),
  output_file = NULL,
  verbose = TRUE
)
```

## Arguments

- se_proc:

  SummarizedExperiment with intensity data

- DEPs_results:

  DataFrame with differential expression results (must have an 'Assay'
  column)

- assay_name:

  Name of the assay to cluster. It is also the value used to filter
  `DEPs_results` by its `Assay` column, so the two must agree. With the
  default `NULL` it is taken from `DEPs_results$Assay`, which records
  the assay the differential abundance was computed on; supply it
  explicitly when those results cover more than one assay.

- filter_mode:

  Filtering mode: "any" (significant in at least one comparison), "all"
  (ALL features, with no significance filtering), "specific"
  (significant in the comparison given by `comparison`)

- alpha:

  Significance threshold (default: 0.05)

- comparison:

  Specific comparison (for filter_mode="specific")

- condition_order:

  Condition order for the profiles

- aggregate:

  Aggregation method: "median" or "mean"

- c_range:

  Range of cluster numbers to evaluate

- auto_select_c:

  Automatic selection of the number of clusters (default: TRUE)

- c:

  Fixed number of clusters (if auto_select_c=FALSE)

- selection_method:

  Selection method: "xb", "consensus", "elbow"

- min_membership:

  Minimum membership threshold for inclusion in the output

- seed:

  Integer. Seed for the Mfuzz clustering itself. Fuzzy c-means starts
  from a random partition, so the cluster labels and the memberships
  depend on it; fix it to make a run reproducible. The caller's RNG
  state is restored on exit. Default 42.

- seeds:

  Numeric vector of seeds used when `auto_select_c = TRUE`: each
  candidate number of clusters is fitted once per seed and the metrics
  are averaged, so that the choice of `c` does not hang on a single
  partition. Default `c(42, 123, 456)`.

- output_file:

  Path of the output parquet file. `NULL` by default, which writes
  nothing to disk; the long-format data is returned regardless in the
  `long_output` element of the result.

- verbose:

  Show progress messages

## Value

List (invisible) with the clustering results: `optimal_c`, `m`,
`conditions`, feature counts, `selection_metrics`, `cluster_counts`, the
Mfuzz `cl` object, the standardised `eset_std` and `long_output`, the
long-format data.frame that is written when `output_file` is given.

## Examples

``` r
if (requireNamespace("Mfuzz", quietly = TRUE) &&
    requireNamespace("e1071", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)

  # auto_select_c = TRUE would cluster once per value of c_range
  pp <- pattern_profiler_analysis(
    res$se_proc,
    DEPs_results    = res$DEPs_results,
    assay_name      = "Impseqrob_min",
    filter_mode     = "any",
    condition_order = c("A", "B", "D"),
    auto_select_c   = FALSE,
    c               = 3,
    verbose         = FALSE
  )
  print(pp$cluster_counts)
  print(head(pp$long_output))
}
#>   hard_assignment Freq
#> 1               1  401
#> 2               2   83
#> 3               3  310
#>   FeatureID Cluster Membership        A        B         D
#> 1    P32610       1   0.977237 0.638391 0.514077 -1.152468
#> 2    Q04182       1   0.977214 0.638649 0.513799 -1.152448
#> 3    P21375       1   0.977183 0.638171 0.514313 -1.152484
#> 4    P36090       1   0.977170 0.638774 0.513665 -1.152439
#> 5    P21147       1   0.977150 0.638096 0.514394 -1.152490
#> 6    Q01476       1   0.977112 0.638027 0.514468 -1.152495
```
