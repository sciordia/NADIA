# Batch Correction with BERT

Applies BERT batch effect correction to a normalized assay in a
SummarizedExperiment. BERT uses hierarchical tree decomposition to
handle missing values, making it suitable for pre-imputation batch
correction on protein-level proteomics data.

## Usage

``` r
batch_correct_proteomics(
  se,
  assay_name,
  batch_column = "Batch",
  corrected_assay_name = "BERT",
  algorithm = "ComBat",
  ComBat_mode = 1,
  covariates = NULL,
  qualitycontrol = FALSE,
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment with a normalized assay.

- assay_name:

  Character. Name of the input assay to correct (e.g., "cycloess").

- batch_column:

  Character. Column in colData(se) containing batch assignments (default
  "Batch").

- corrected_assay_name:

  Character. Name for the new corrected assay added to the SE (default
  "BERT").

- algorithm:

  Character: "ComBat" (default), "limma", or "ref".

- ComBat_mode:

  Integer 1-4 controlling ComBat behavior: 1 = parametric +
  mean+variance, 2 = parametric + mean-only, 3 = non-parametric +
  mean+variance, 4 = non-parametric + mean-only.

- covariates:

  Character vector of column names from colData(se) to protect during
  batch correction (default NULL). These are mapped to BERT's Cov_1,
  Cov_2, ... format internally. BERT uses those columns directly as the
  design matrix and encodes nothing itself, so categorical columns are
  converted to indicator variables here; numeric columns are passed
  through unchanged and act as continuous covariates. IMPORTANT: for
  ComBat/limma, include here the biological variable of interest (e.g.
  the condition) in order to PRESERVE it; otherwise ComBat removes all
  batch-associated variance and can erase biological signal when
  condition and batch are confounded.

- qualitycontrol:

  Logical. Compute ASW (Average Silhouette Width) quality metrics for
  raw vs corrected data (default FALSE).

- verbose:

  Logical (default TRUE).

## Value

SummarizedExperiment with the new corrected assay added. If BERT returns
fewer features than the input, the SE is subsetted to match and a
warning is issued.

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

# BERT removes the batch effect and adds a "BERT" assay, leaving the
# original one untouched
se <- batch_correct_proteomics(res$se_proc, assay_name = "Impseqrob_min",
                               batch_column = "Batch", verbose = FALSE)
#> Warning: batch_correct_proteomics: 'covariates = NULL' with algorithm='ComBat'. ComBat/limma remove ALL batch-associated variance; if the biological condition is (partially) confounded with the batch, real biological signal will be lost. Passing the condition variable in 'covariates' (batch_covariates) is recommended in order to preserve it.
SummarizedExperiment::assayNames(se)
#> [1] "raw"           "log2"          "cycloess"      "Impseqrob_min"
#> [5] "BERT"         

# On a real experiment you would normally pass the biological variable of
# interest in `covariates`, so that ComBat preserves it instead of removing
# it along with the batch effect. `covariates = "Condition"` works on this
# example: the assay is already imputed, so each batch-by-condition cell has
# its two values. Before imputation, with only two samples per cell, BERT
# needs two non-missing values per cell and any protein with a gap there
# stops being adjustable.
```
