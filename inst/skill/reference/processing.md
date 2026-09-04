# Processing: `process_proteomics()` and its stages

One coordinator runs four stages in order:

**Normalisation → (optional) batch correction → imputation → differential
abundance.**

Pattern Profiler is **not** part of it; it is a separate step run afterwards
(see `reference/visualisation.md`).

## The arguments that matter

```r
res <- process_proteomics(
  preprocessing,                    # a proteomics_data object

  # --- filtering ---
  min_reps_filter   = NULL,         # min replicates with a value, per group
  min_groups_filter = 1,            # min groups meeting that

  # --- normalisation ---
  norm_method = "cycloess",

  # --- batch correction (off by default) ---
  batch_correct     = FALSE,
  batch_column      = "Batch",
  batch_algorithm   = "ComBat",     # or "limma"
  batch_covariates  = NULL,         # PASS THE CONDITION HERE. See below.
  covariate_df      = NULL,         # data frame carrying the batch column

  # --- imputation ---
  imp_method  = "combo",
  mar_method  = "Impseqrob",
  mnar_method = "min",
  max_na_prop = NULL,               # NULL = no filtering by missingness
  method_args = list(),             # per-method extra arguments

  # --- differential abundance ---
  de_method       = "limma",        # or "limpa"
  comparisons     = NULL,           # NULL = all pairwise
  control         = NULL,           # or: everything against this condition
  alpha           = 0.05,
  logFC_threshold = 0,
  covariate_column     = NULL,      # adjust the model for a covariate
  bio_replicate_column = NULL,      # blocking factor (paired designs)

  # --- export (see reference/export.md) ---
  export_dir = NULL,                # NULL = nothing written to disk

  verbose = TRUE
)
```

### Filtering

`min_reps_filter` / `min_groups_filter` drop proteins not seen often enough.
Applied at normalisation, before anything is imputed. Leave `min_reps_filter`
at `NULL` unless you have a reason; be aware that filtering hard here removes
exactly the presence/absence proteins that DIA experiments are often about.

### Choosing the contrasts

- `comparisons = NULL` → every pairwise contrast. Three conditions give three.
- `comparisons = c("B-A", "D-A")` → only those.
- `control = "A"` → everything against A.

Format is always `"numerator-denominator"`; a positive `logFC` is higher in the
numerator.

### Batch correction

Off by default. Turn it on only when you have a real batch structure and
evidence that it matters — check with a PCA coloured by batch, or the PVCA
decomposition:

```r
pca_covariates_plot(res$se_proc, assay_name = "cycloess",
                    covariates = c("Condition", "Batch"))
```

Then:

```r
res <- process_proteomics(
  prep,
  covariate_df     = my_batch_table,      # must contain the batch column
  batch_correct    = TRUE,
  batch_column     = "Batch",
  batch_algorithm  = "ComBat",
  batch_covariates = "Condition"          # <- protects the biology
)
```

**`batch_covariates` is not optional in practice.** With `NULL`, ComBat removes
every batch-associated component, and when condition and batch are partly
confounded that includes the biological signal. NADIA emits a warning; treat it
as an error unless you can justify otherwise.

Batch correction adds a `"BERT"` assay and imputation then runs on it instead
of on the normalised assay. Features that cannot be fitted (fewer than two
finite observations, or zero within-batch variance) are set aside and
reinserted unadjusted rather than dropped.

To run it on its own: `batch_correct_proteomics(se, assay_name, batch_column,
algorithm, covariates, ...)`.

### limma or limpa

`de_method = "limma"` is the default and works with every imputation method.

`de_method = "limpa"` must be paired with `imp_method = "limpa"`: that method
attaches an `EList` with precision weights which the DE step then uses. Using
one without the other wastes the point of it.

`eBayes_trend` and `eBayes_robust` default to `NULL` and are resolved from
`de_method` (TRUE for limma, FALSE for limpa). Override only deliberately.

## Running the stages separately

The coordinator is a convenience, not a requirement:

```r
norm <- normalize_proteomics(prep, norm_method = "cycloess")
se   <- norm$se
se   <- batch_correct_proteomics(se, assay_name = "cycloess", batch_column = "Batch",
                                 covariates = "Condition")
se   <- impute_proteomics(se, assay_name = "BERT", imp_method = "combo")
de   <- de_analysis_proteomics(se, assay_name = "Impseqrob_min", de_method = "limma")
```

Useful when you want to try several imputations on one normalised object
without repeating the normalisation.

## Checking the output before believing it

```r
SummarizedExperiment::assayNames(res$se_proc)   # which stages ran
dim(res$se_proc)                                # proteins x samples after filtering
table(res$DEPs_results$Comparison, res$DEPs_results$Change)

# How much of each comparison rests on imputed values?
with(res$DEPs_results[res$DEPs_results$Change != "No Change", ],
     summary(MissComp))
```

A significant set whose `MissComp` is mostly high is an imputation artefact
story, not a biology story. Report the distribution, not just the count.

`S4Vectors::metadata(res$se_proc)` also keeps `mar_mask` and `mnar_mask`: which
cells were filled, and by which branch.
