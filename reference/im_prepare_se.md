# Build a SummarizedExperiment for imputation benchmarking

Convenience function that creates a SE from a `proteomics_data` object
and applies the winning normalization method (from normalization
benchmarking). The result is ready for
[`imputation_metrics()`](https://sciordia.github.io/NADIA/reference/imputation_metrics.md).

## Usage

``` r
im_prepare_se(
  preprocessing,
  norm_method = NULL,
  pc1_rank = NULL,
  min_reps = NULL,
  min_groups = 1,
  covariate_df = NULL,
  norm_method_args = list(),
  include_log2 = TRUE,
  verbose = TRUE
)
```

## Arguments

- preprocessing:

  A `proteomics_data` object, from any of the `preprocess_*()`
  functions.

- norm_method:

  Character scalar. Normalization method name (e.g. `"cycloess"`). If
  provided, used directly. Default `NULL`.

- pc1_rank:

  data.frame from
  [`nm_rank_pc1()`](https://sciordia.github.io/NADIA/reference/nm_rank_pc1.md)
  or `normalization_metrics()$pc1_rank`. The first row's `Method` column
  is used as winner. Ignored if `norm_method` is provided. Default
  `NULL`.

- min_reps:

  Minimum replicates with non-NA values per group for protein filtering.
  If NULL, auto-computed as half the smallest group. Default `NULL`.

- min_groups:

  Minimum groups meeting `min_reps` (default: 1).

- covariate_df:

  Optional covariate data.frame for paired designs (must contain a
  `Column` column). Default `NULL`.

- norm_method_args:

  Named list of per-method arguments for normalization. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- include_log2:

  Logical. Include the baseline `"log2"` assay in the returned SE.
  Default `TRUE`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

SummarizedExperiment with assays `"log2"` (optional) + winner method.

## Details

Internally calls
[`nm_prepare_se()`](https://sciordia.github.io/NADIA/reference/nm_prepare_se.md)
to build the baseline SE with assay `"log2"`, then applies the chosen
normalization via
[`.nm_dispatch_normalization()`](https://sciordia.github.io/NADIA/reference/dot-nm_dispatch_normalization.md).

## Examples

``` r
data(nadia_dia)

# Option A: name the normalization to benchmark imputation on
se <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)
SummarizedExperiment::assayNames(se)
#> [1] "log2"     "cycloess"

# Option B: let the winner of normalization_metrics() decide. That function
# takes a SummarizedExperiment, so the baseline is built first.
se_base <- nm_prepare_se(nadia_dia, verbose = FALSE)
nm <- normalization_metrics(se_base,
                            methods = c("cycloess", "quantile", "MAD"),
                            plots = "pc1_ranking", verbose = FALSE)
nm$pc1_rank
#>     Method PC1_VarPct Rank
#> 1 cycloess   72.57836    1
#> 2 quantile   65.48778    2
#> 3      MAD   64.86712    3
#> 4     log2   56.03859    4
se2 <- im_prepare_se(nadia_dia, pc1_rank = nm$pc1_rank, verbose = FALSE)
SummarizedExperiment::assayNames(se2)
#> [1] "log2"     "cycloess"
```
