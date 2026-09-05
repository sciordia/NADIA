# Process proteomics data

Complete processing pipeline: filtering, normalization, imputation and
differential expression analysis. Coordinates sub-modules
Normalization.R, Imputation.R, and DEAnalysis.R.

## Usage

``` r
process_proteomics(
  preprocessing,
  export_dir = NULL,
  min_reps_filter = NULL,
  min_groups_filter = 1,
  norm_method = "cycloess",
  cyclic_loess_method = "fast",
  cyclic_loess_iterations = 3,
  cyclic_loess_span = 0.7,
  batch_correct = FALSE,
  batch_column = "Batch",
  batch_algorithm = "ComBat",
  batch_ComBat_mode = 1,
  batch_covariates = NULL,
  batch_qualitycontrol = FALSE,
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
  comparisons = NULL,
  control = NULL,
  logFC_threshold = 0,
  alpha = 0.05,
  eBayes_trend = NULL,
  eBayes_robust = NULL,
  de_method = "limma",
  covariate_df = NULL,
  covariate_column = NULL,
  bio_replicate_column = NULL,
  export_normalized = TRUE,
  export_imputed = TRUE,
  export_format = "tsv",
  export_volcano = TRUE,
  export_boxplot = TRUE,
  export_pca = TRUE,
  verbose = TRUE
)
```

## Arguments

- preprocessing:

  A `proteomics_data` object, as returned by
  [`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md),
  [`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md),
  [`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md)
  or
  [`preprocess_lfq()`](https://sciordia.github.io/NADIA/reference/preprocess_lfq.md).

- export_dir:

  Output directory for exported files. Defaults to `NULL`, which writes
  nothing to disk; pass a path to enable the exports controlled by
  `export_normalized`, `export_imputed`, `export_volcano`,
  `export_boxplot` and `export_pca`.

- min_reps_filter:

  Minimum replicates for filtering. If NULL, auto-computed

- min_groups_filter:

  Minimum groups for filtering (default: 1)

- norm_method:

  Normalization method passed to normalize_proteomics() (default:
  "cycloess"). See normalize_proteomics() for all 13 options.

- cyclic_loess_method:

  Cyclic Loess method: "fast" or "pairs" (default: "fast")

- cyclic_loess_iterations:

  Number of iterations for Cyclic Loess (default: 3)

- cyclic_loess_span:

  Span parameter for Cyclic Loess (default: 0.7)

- batch_correct:

  Logical. Apply BERT batch correction after normalization (default:
  FALSE). Requires covariate_df with batch_column.

- batch_column:

  Column in covariate_df containing batch assignments (default:
  "Batch"). Must have \>= 2 unique values.

- batch_algorithm:

  Batch correction algorithm: "ComBat" (default), "limma", or "ref"

- batch_ComBat_mode:

  Integer 1-4 for ComBat parametric/mean-only settings (default: 1)

- batch_covariates:

  Character vector of colData column names to protect during BERT batch
  correction (default: NULL). Categorical columns are encoded as
  indicator variables automatically; numeric columns are passed through
  unchanged and therefore act as continuous covariates.

- batch_qualitycontrol:

  Logical. Compute ASW quality metrics (default: FALSE)

- imp_method:

  Imputation method (default: "combo"). See impute_proteomics() for all
  options.

- mar_method:

  MAR method for combo mode (default: "Impseqrob")

- mnar_method:

  MNAR method for combo mode (default: "min")

- prop_na_mnar:

  NA proportion for MNAR classification (default: 0.51)

- prop_present_mar:

  Present proportion for MAR (default: 0.5)

- min_present_mar:

  Minimum present values for MAR (default: 1)

- require_n_conditions:

  Required conditions with presence (default: 1)

- max_na_prop:

  Maximum NA proportion above which a protein is removed before
  imputation. `NULL` (the default) disables the filter.

- method_args:

  Named list of per-method argument lists

- with_value:

  Constant value for imp_method="with"

- comparisons:

  Comparisons for DE. If NULL, generates all pairwise

- control:

  Control condition. If NULL, compares all

- logFC_threshold:

  LogFC threshold for significance (default: 0)

- alpha:

  Adjusted p-value threshold (default: 0.05)

- eBayes_trend:

  Use trend estimation in eBayes. If NULL (default), it is resolved from
  de_method: TRUE for "limma", FALSE for "limpa".

- eBayes_robust:

  Use robust estimation in eBayes. If NULL (default), it is resolved
  from de_method: TRUE for "limma", FALSE for "limpa".

- de_method:

  DE method: "limma" (default) or "limpa" (probabilistic, requires
  imp_method="limpa")

- covariate_df:

  Data frame with Column + covariate column(s) for paired/blocked design
  (default: NULL)

- covariate_column:

  Name(s) of the covariate column(s) for the DE model. Single string
  (e.g., "Subject") or character vector (e.g., c("Subject", "Batch")).
  Default: NULL

- bio_replicate_column:

  Column name in covariate_df identifying biological replicates (e.g.,
  "Patient", "Subject"). Used with limma::duplicateCorrelation() to
  account for technical replicates or paired designs via random effect
  blocking. Default: NULL

- export_normalized:

  Export normalized matrix (default: TRUE)

- export_imputed:

  Export imputed matrix (default: TRUE)

- export_format:

  Export format: "tsv", "parquet", or "both" (default: "tsv")

- export_volcano:

  Export VolcanoPlot_Input (default: TRUE)

- export_boxplot:

  Export BoxPlot_Input (default: TRUE)

- export_pca:

  Export PCA_Input (default: TRUE)

- verbose:

  Print progress messages (default: TRUE)

## Value

List of class proteomics_result with:

- se_proc: Processed SummarizedExperiment with all assays

- DEPs_results: Data frame with differential expression results

- BoxPlot_Input: Long-format data frame with one row per protein, sample
  and assay (columns Column, Assay, Intensity, Condition, Replicate,
  Protein.IDs). Feed it to
  [`boxplot_highchart_list()`](https://sciordia.github.io/NADIA/reference/boxplot_highchart_list.md).

- PCA_Input: Long-format data frame with the imputed intensities plus
  one adjP\_\* column per comparison and a sig_any flag. Feed it to
  [`pca_highchart_list()`](https://sciordia.github.io/NADIA/reference/pca_highchart_list.md),
  [`proteomics_heatmap()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap.md)
  or
  [`proteomics_heatmap_list()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap_list.md).

- comparisons: Comparisons performed

- parameters: Parameters used

`BoxPlot_Input` and `PCA_Input` are always returned, whether or not
`export_dir` is given: the plotting layer needs them, and requiring a
round trip through disk to obtain them would be gratuitous.

## Details

The input is a `proteomics_data` object, which every preprocessing
function returns with the same three elements, so the pipeline is
indifferent to the acquisition and the search engine the data came from.

## Examples

``` r
data(nadia_dia)

# Normalization, imputation and differential expression in one call
res <- process_proteomics(nadia_dia, verbose = FALSE)
res
#> === Proteomics Processing Result ===
#> 
#> SummarizedExperiment:
#>   - Proteins: 1997 
#>   - Samples: 12 
#>   - Assays: raw, log2, cycloess, Impseqrob_min 
#>   - Conditions: A, B, D 
#> 
#> Differential Results:
#>   - Total rows: 5991 
#>   - Comparisons: B-A, D-A, D-B 
#>     B-A: Up=325, Down=172
#>     D-A: Up=374, Down=383
#>     D-B: Up=286, Down=382
#> 
#> Parameters:
#>   - Normalization: cycloess 
#>     - Cyclic Loess method: fast 
#>     - Cyclic Loess iterations: 3 
#>     - Cyclic Loess span: 0.7 
#>   - Imputation: combo 
#>     - MAR method: Impseqrob 
#>     - MNAR method: min 
#>   - Alpha: 0.05 
#>   - logFC threshold: 0 
#>   - Output directory: (no export) 

SummarizedExperiment::assayNames(res$se_proc)
#> [1] "raw"           "log2"          "cycloess"      "Impseqrob_min"
table(res$DEPs_results$Comparison, res$DEPs_results$Change)
#>      
#>         Up Down No Change
#>   B-A  325  172      1500
#>   D-A  374  383      1240
#>   D-B  286  382      1329

# The plotting layer is fed straight from the result, no files involved
names(res$BoxPlot_Input)
#> [1] "Column"      "Assay"       "Intensity"   "Condition"   "Replicate"  
#> [6] "Protein.IDs"
names(res$PCA_Input)
#> [1] "SampleID"  "FeatureID" "Intensity" "Condition" "Replicate" "adjP_B-A" 
#> [7] "adjP_D-A"  "adjP_D-B"  "sig_any"  

# Another combination: quantile normalization and a single imputation method
res2 <- process_proteomics(nadia_dia, norm_method = "quantile",
                           imp_method = "QRILC", verbose = FALSE)
SummarizedExperiment::assayNames(res2$se_proc)
#> [1] "raw"      "log2"     "quantile" "QRILC"   
```
