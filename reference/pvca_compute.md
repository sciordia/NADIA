# Compute PVCA variance components

Performs Principal Variance Component Analysis on a SummarizedExperiment
object. Decomposes total variance into contributions from each specified
factor, their pairwise interactions, and residual variance.

## Usage

``` r
pvca_compute(
  se,
  assay_name = NULL,
  factors,
  pca_threshold = 0.6,
  variance_threshold = 0.01,
  na_action = "complete",
  fill_value = -1,
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment object.

- assay_name:

  Character. Name of the assay to use. If NULL, uses the second assay
  (typically the normalized one).

- factors:

  Character vector. Column names in colData(se) to include as variance
  components. All specified factors must exist in colData.

- pca_threshold:

  Numeric (0-1). Minimum cumulative proportion of variance explained by
  retained principal components (default 0.6).

- variance_threshold:

  Numeric (0-1). Factors explaining less than this proportion are
  grouped into a "Below X%" category (default 0.01).

- na_action:

  Character. How to handle NAs: "complete" (default) removes rows with
  any NA, "fill" replaces NAs with fill_value, "none" assumes no NAs
  (errors if present).

- fill_value:

  Numeric. Replacement for NAs when na_action = "fill" (default -1).

- verbose:

  Logical. Print progress messages (default TRUE).

## Value

data.frame with columns:

- label:

  Factor name, interaction term, "Below X%", or "resid"

- weights:

  Weighted average proportion of variance (0-1)

## Details

Uses lme4 directly to fit mixed models, which allows automatic detection
and exclusion of interaction terms where the number of levels equals or
exceeds the number of observations (e.g., Condition:Patient in paired
designs with one observation per patient-condition combination).

## Examples

``` r
data(nadia_dia)

# nadia_dia has no batch information, so the example builds a plausible one:
# two digestion batches crossed with the three conditions.
cov <- data.frame(Column = nadia_dia$metadata$Coding,
                  Batch = rep(c("b1", "b2"),
                              length.out = nrow(nadia_dia$metadata)),
                  stringsAsFactors = FALSE)
res <- process_proteomics(nadia_dia, covariate_df = cov, verbose = FALSE)

# Weight of each factor, and of their interaction, on the total variance
vc <- pvca_compute(res$se_proc, assay_name = "Impseqrob_min",
                   factors = c("Condition", "Batch"), verbose = FALSE)
#> boundary (singular) fit: see help('isSingular')
vc
#>       label      weights
#> 1 Condition 9.999576e-01
#> 2  Below 1% 4.241289e-05
```
