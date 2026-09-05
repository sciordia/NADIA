# Impute proteomics data

Complete imputation pipeline with 20 methods. Four pathways:

## Usage

``` r
impute_proteomics(
  se,
  normalized_assay_name = "cycloess",
  imputed_assay_name = NULL,
  imp_method = "combo",
  mar_method = "Impseqrob",
  mnar_method = "min",
  prop_na_mnar = 0.51,
  prop_present_mar = 0.5,
  min_present_mar = 1,
  require_n_conditions = 1,
  max_na_prop = NULL,
  method_args = list(),
  with_value = NA_real_,
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment with normalized assay

- normalized_assay_name:

  Name of the normalized assay to use (default: "cycloess")

- imputed_assay_name:

  Name for the imputed assay. If NULL, it is built from the method:
  `combo` gives `<mar_method>_<mnar_method>`, `softHybrid` gives
  `softHybrid_<mar_method>_<mnar_method>`, and any single method gives
  its own name.

- imp_method:

  Imputation method (default: "combo"). One of: "combo", "softHybrid",
  "bpca", "knn", "mice", "missForest", "Impseq", "Impseqrob", "QRILC",
  "MLE", "MinDet", "MinProb", "PI", "min", "halfmin", "zero", "nbavg",
  "with", "limpa", "none". `"halfmin"` is the DIA-NN half-minimum
  recipe: since the assay is on the log2 scale, halving the intensity
  means one unit below the global observed minimum.

- mar_method:

  MAR method for combo/softHybrid mode (default: "Impseqrob")

- mnar_method:

  MNAR method for combo/softHybrid mode (default: "min")

- prop_na_mnar:

  NA proportion threshold for MNAR classification (default: 0.51)

- prop_present_mar:

  Present proportion for MAR (default: 0.5)

- min_present_mar:

  Minimum present values for MAR (default: 1)

- require_n_conditions:

  Number of conditions required with presence (default: 1)

- max_na_prop:

  Maximum NA proportion above which a protein is removed before
  imputation. `NULL` (the default) disables the filter, so every method
  receives the same proteins.

- method_args:

  Named list of per-method argument lists (e.g., list(knn = list(k =
  10), bpca = list(nPcs = 3)))

- with_value:

  Constant value for imp_method="with" (default: NA_real\_)

- verbose:

  Print progress messages (default: TRUE)

## Value

List with:

- se: SummarizedExperiment with imputed assay added

- x_imputed: Raw imputed matrix (for direct export)

- prefilter_summary: Pre-filtering summary

- imputation_summary: Imputation statistics

- mnar_mask: MNAR mask matrix (combo only, NULL otherwise)

- mar_mask: MAR mask matrix (combo only, NULL otherwise)

## Details

- `imp_method = "none"`: no imputation

- `imp_method = "combo"`: two-stage MAR+MNAR with binary classification

- `imp_method = "softHybrid"`: sigmoid-weighted MAR+MNAR blend
  (continuous)

- Any other method: apply to all NAs (no MAR/MNAR distinction)

## Examples

``` r
data(nadia_dia)
norm <- normalize_proteomics(nadia_dia, norm_method = "cycloess",
                             verbose = FALSE)

# Default: two-stage combo, Impseqrob for MAR and min for MNAR
imp <- impute_proteomics(norm$se, normalized_assay_name = "cycloess",
                         verbose = FALSE)
SummarizedExperiment::assayNames(imp$se)
#> [1] "raw"           "log2"          "cycloess"      "Impseqrob_min"
sum(is.na(SummarizedExperiment::assay(imp$se, "Impseqrob_min")))
#> [1] 0

# Any single method works too, and per-method options go in method_args
imp2 <- impute_proteomics(norm$se, normalized_assay_name = "cycloess",
                          imp_method = "knn",
                          method_args = list(knn = list(k = 15)),
                          verbose = FALSE)
SummarizedExperiment::assayNames(imp2$se)
#> [1] "raw"      "log2"     "cycloess" "knn"     
```
