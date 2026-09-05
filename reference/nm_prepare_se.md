# Prepare a SummarizedExperiment from a proteomics_data object

Convenience wrapper that extracts metadata and protein data from a
`proteomics_data` object (from any of the preprocess\_\*() functions),
performs zero-to-NA conversion, protein filtering, and returns a SE with
assays `"raw"` and `"log2"` – ready for
[`nm_run_normalizations()`](https://sciordia.github.io/NADIA/reference/nm_run_normalizations.md)
or `normalization_metrics(..., methods = "all")`.

## Usage

``` r
nm_prepare_se(
  preprocessing,
  min_reps = NULL,
  min_groups = 1,
  covariate_df = NULL,
  verbose = TRUE
)
```

## Arguments

- preprocessing:

  A `proteomics_data` object, from any of the `preprocess_*()`
  functions.

- min_reps:

  Minimum replicates with non-NA values per group for protein filtering.
  If NULL, auto-computed as half the smallest group. Default `NULL`.

- min_groups:

  Minimum groups meeting `min_reps` (default: 1).

- covariate_df:

  Optional covariate data.frame for paired designs (must contain a
  `Column` column). Default `NULL`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

SummarizedExperiment with assays `"raw"` and `"log2"`.

## Details

Internally calls
[`normalize_proteomics()`](https://sciordia.github.io/NADIA/reference/normalize_proteomics.md)
with `norm_method = "log2"` (no additional normalization).

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
dim(se)
#> [1] 1997   12
SummarizedExperiment::assayNames(se)
#> [1] "raw"  "log2"
```
