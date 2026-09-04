# Choosing the normalisation and imputation methods

This is the stage that changes the answer. On the bundled example, four
combinations of the same data, the same model and the same thresholds give:

| Pipeline | Up | Down | No Change |
|---|---|---|---|
| `cycloess` + `combo` | 985 | 937 | 4069 |
| `cycloess` + `knn` | 872 | 179 | 4940 |
| `quantile` + `combo` | 902 | 2825 | 2264 |
| `log2` + `min` | 761 | 778 | 4452 |

Down ranges from 179 to 2,825 — a sixteen-fold difference produced entirely by
method choice. **Never present one pipeline's output without saying which
pipeline it was, and never pick one because it gives more hits.**

## Normalisation: 13 methods

```
log2  log2Norm  GlobalMedian  GlobalMean  eqmedians  medianNorm  meanNorm
quantile  quantile.robust  Rlr  MAD  cycloess  vsn
```

Default `cycloess`. Notes that matter:

- `log2` adds no new assay: it log-transforms and stops. Use it as the "no
  normalisation" control.
- `GlobalMedian` and `GlobalMean` normalise by **column sum**, not by column
  median/mean, and differ only by a global additive constant. Downstream they
  are effectively the same method; do not report both as independent evidence.
- `vsn`, `Rlr` need `Suggests` packages (`vsn`, `MASS`).
- `quantile` uses `limma::normalizeQuantiles` for consistent NA handling.

## Imputation: 20 methods, three kinds

**Hybrid (2), the interesting ones:**

- `combo` (default) — two stages. `mar_method` fills values assumed missing at
  random, `mnar_method` fills the rest. Defaults `Impseqrob` + `min`. MNAR
  parameters are estimated on the *original* matrix, not the MAR-filled one, so
  the low-end estimate is not biased upward.
- `softHybrid` — a sigmoid-weighted blend of a MAR and an MNAR method rather
  than a hard split.

**MAR stage candidates:** `bpca knn mice missForest Impseq Impseqrob MLE none`

**MNAR stage candidates:** `QRILC MinDet MinProb PI min halfmin zero with none`

**Single methods**, applied to every missing value alike: any of the above plus
`nbavg`, `limpa`, and `none`.

No extra dependency: `none zero min halfmin MinDet PI nbavg with`.
Needs a `Suggests` package: `bpca knn mice missForest Impseq Impseqrob QRILC
MLE MinProb limpa`.

`limpa` is different in kind — a model-based method that also carries precision
weights through to the DE step. Use it with `de_method = "limpa"`.

## First: is there anything to impute?

Before comparing imputation methods, count the missing cells. If there are
none, or a handful, the whole imputation question is moot — every method
returns the input — and `imputation_metrics()` would rank methods on an
artificial problem. Say so and move on to normalisation, which still matters.
The TMT example shipped with the package is exactly this case: one missing cell
in 48,000.

## How to choose, in order of preference

**1. If the experiment has a known ground truth (a spike-in), benchmark it.**
That is the only way to measure which pipeline recovers the real changes. See
`reference/validation.md`, `benchmarking_proteomics()` and
`benchmarking_multiple()`.

**2. Otherwise, rank on data-derived metrics.** Two orchestrators, no ground
truth needed:

```r
se <- normalize_proteomics(prep, norm_method = "cycloess")$se   # or res$se_proc

nm <- normalization_metrics(
  se,
  methods = c("log2Norm", "GlobalMedian", "GlobalMean", "eqmedians", "vsn",
              "medianNorm", "meanNorm", "quantile", "Rlr", "MAD", "cycloess",
              "quantile.robust"),          # required: see reference/validation.md
  output_dir = "results/norm")
nm$final_rank        # rank 1 = best

im <- imputation_metrics(se, assay_name = "cycloess", output_dir = "results/imp")
im$metrics_table     # Rank_Mean over NRMSE, SOR, PSS, ACC_OI
```

`normalization_metrics()` compares those 12 methods on within-group variability
(PCV/PMAD/PEV), intragroup correlation and group separation. **Omit `methods`
and it ranks the assays already in the object instead**, which is not a method
comparison at all — `reference/validation.md` explains why that matters.
`imputation_metrics()` simulates missingness on complete rows and measures
recovery for 16 methods.

**3. Only then fall back to the defaults**, and say that is what you did.

## Reading the missingness first

Choose the *kind* of method by looking at the data, before running anything:

```r
se <- res$se_proc
x  <- SummarizedExperiment::assay(se, "cycloess")     # normalised, pre-imputation
mean(is.na(x))                                        # overall missingness

# Is missingness concentrated in the low-intensity range? -> MNAR component
obs <- rowMeans(x, na.rm = TRUE)
plot(obs, rowMeans(is.na(x)))
```

A clear negative relationship between mean intensity and missing rate means a
real MNAR component, and a purely MAR method (`knn`, `missForest`) will pull
those values up towards the observed mean and erase the effect. That is exactly
the `cycloess` + `knn` row in the table above: 179 Down against 937.

If missingness looks flat across intensity, the MNAR stage matters less and a
MAR method alone is defensible.

## What NADIA leaves alone

There is no automatic method selection. `process_proteomics()` runs what you
ask for. Choosing is the analyst's job, and this skill's job is to make that
choice explicit and justified — not to hide it behind a default.
