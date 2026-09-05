# Normalize proteomics data

Complete normalization pipeline: zero-to-NA conversion, protein
filtering by group presence, SummarizedExperiment creation, and
normalization with the selected method. Output is always in log2 scale.

## Usage

``` r
normalize_proteomics(
  data,
  metadata = NULL,
  min_reps = NULL,
  min_groups = 1,
  norm_method = "cycloess",
  cyclic_loess_method = c("fast", "pairs"),
  cyclic_loess_iterations = 3,
  cyclic_loess_span = 0.7,
  verbose = TRUE
)
```

## Arguments

- data:

  Either a `proteomics_data` object – the output of
  [`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md),
  [`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md),
  [`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md)
  or
  [`preprocess_lfq()`](https://sciordia.github.io/NADIA/reference/preprocess_lfq.md)
  – in which case `metadata` is derived from it and must be left `NULL`,
  or a data frame with ProteinGroups, GeneNames, UniqPepts and the
  intensity columns.

- metadata:

  Data frame with Column, Condition, Replicate. Only needed when `data`
  is a data frame.

- min_reps:

  Minimum replicates for filtering. NULL = auto: floor(min_group_size /
  2)

- min_groups:

  Minimum groups meeting min_reps (default: 1)

- norm_method:

  Normalization method (default: "cycloess"). One of:

  - Group A (input: raw intensities): "log2Norm", "GlobalMedian",
    "GlobalMean", "eqmedians", "vsn"

  - Group B (input: log2 assay): "log2" (no extra normalization),
    "quantile", "Rlr", "MAD", "cycloess", "medianNorm", "meanNorm",
    "quantile.robust"

- cyclic_loess_method:

  Cyclic Loess method: "fast" or "pairs" (default: "fast")

- cyclic_loess_iterations:

  Number of iterations for Cyclic Loess (default: 3)

- cyclic_loess_span:

  Span parameter for Cyclic Loess (default: 0.7)

- verbose:

  Print progress messages (default: TRUE)

## Value

List with:

- se: SummarizedExperiment with assays raw, log2, and \<norm_method\>
  (assays raw + log2 only when norm_method = "log2")

- filter_summary: Filtering summary

- na_overview: NA statistics

## Examples

``` r
data(nadia_dia)

# The object returned by preprocess_*() is accepted directly
norm <- normalize_proteomics(nadia_dia, norm_method = "cycloess",
                             verbose = FALSE)
SummarizedExperiment::assayNames(norm$se)
#> [1] "raw"      "log2"     "cycloess"

# Zeros become NA and low-coverage proteins are filtered out
norm$filter_summary
#> $n_total
#> [1] 2000
#> 
#> $n_keep
#> [1] 1997
#> 
#> $n_drop
#> [1] 3
#> 
#> $min_reps
#> [1] 2
#> 
#> $min_groups
#> [1] 1
#> 

# Any of the 13 methods can be used
norm2 <- normalize_proteomics(nadia_dia, norm_method = "quantile",
                              verbose = FALSE)
SummarizedExperiment::assayNames(norm2$se)
#> [1] "raw"      "log2"     "quantile"
```
