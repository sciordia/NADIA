# Perform differential expression analysis on proteomics data

Coordinates comparison specification, limma analysis, and result
extraction.

## Usage

``` r
de_analysis_proteomics(
  se,
  assay_name = NULL,
  comparisons = NULL,
  control = NULL,
  logFC_threshold = 0,
  alpha = 0.05,
  p_adj = TRUE,
  eBayes_trend = NULL,
  eBayes_robust = NULL,
  de_method = "limma",
  covariate_column = NULL,
  bio_replicate_column = NULL,
  condition_column = "Condition",
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment with imputed assay (output of impute_proteomics)

- assay_name:

  Assay name to use. If NULL, uses the last available assay

- comparisons:

  Comparisons to perform. If NULL, auto-generates (pairwise or vs
  control)

- control:

  Control condition. If NULL, generates all pairwise comparisons

- logFC_threshold:

  LogFC threshold for significance (default: 0)

- alpha:

  Adjusted p-value threshold (default: 0.05)

- p_adj:

  Use adjusted p-value (default: TRUE)

- eBayes_trend:

  Use trend estimation in eBayes (default: TRUE, recommended for
  proteomics)

- eBayes_robust:

  Use robust estimation in eBayes (default: TRUE, recommended for
  proteomics)

- de_method:

  DE method: "limma" (default) or "limpa" (probabilistic, requires
  imp_method="limpa")

- covariate_column:

  Column name(s) in colData for paired/blocked design. Single string
  (e.g., "Subject") or character vector (e.g., c("Subject", "Batch")).
  Default: NULL

- bio_replicate_column:

  Column name in colData identifying biological replicates (e.g.,
  "Patient", "Subject"). Used with limma::duplicateCorrelation() to
  account for technical replicates or paired designs via random effect
  blocking. Default: NULL

- condition_column:

  Condition column name (default: "Condition")

- verbose:

  Print progress messages (default: TRUE)

## Value

List with:

- DEPs_results: Data frame with differential expression results

- comparisons: Comparisons performed

- assay_name: Assay name used

## Examples

``` r
data(nadia_dia)
norm <- normalize_proteomics(nadia_dia, norm_method = "cycloess",
                             verbose = FALSE)
imp  <- impute_proteomics(norm$se, normalized_assay_name = "cycloess",
                          verbose = FALSE)

# All pairwise comparisons
de <- de_analysis_proteomics(imp$se, assay_name = "Impseqrob_min",
                             verbose = FALSE)
head(de$DEPs_results[, c("Protein.IDs", "Comparison", "logFC", "adj.P.Val")])
#>         Protein.IDs Comparison       logFC    adj.P.Val
#> 1        A0A024RBG1        B-A  0.08293965 7.389228e-01
#> 2        A0A024RBG1        D-A -0.01391145 9.511693e-01
#> 3        A0A024RBG1        D-B -0.09685110 7.151013e-01
#> 4 A0A140T897;P02769        B-A  0.13439560 3.118509e-03
#> 5 A0A140T897;P02769        D-B  0.20481573 7.760613e-05
#> 6 A0A140T897;P02769        D-A  0.33921133 5.778484e-07

# Or every condition against a single control
de2 <- de_analysis_proteomics(imp$se, assay_name = "Impseqrob_min",
                              control = "A", verbose = FALSE)
unique(de2$DEPs_results$Comparison)
#> [1] "B-A" "D-A"
```
