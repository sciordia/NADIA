---
name: nadia
description: Run a differential abundance analysis of bottom-up proteomics data with the R package NADIA — import Spectronaut, DIA-NN or Proteome Discoverer (TMT/LFQ) reports, normalise, correct batch effects, impute missing values, test with limma or limpa, benchmark the method choice, visualise, and export or archive the result. Use whenever the task involves proteomics quantification reports, protein-level differential abundance, DIA missing values, or the NADIA package itself.
---

# Analysing proteomics data with NADIA

NADIA takes a protein quantification report and produces differential abundance
results, treating missing values as information rather than as noise to be
removed. This skill covers the whole pipeline and, more importantly, the
decisions inside it.

## Before anything else

NADIA is an R package. It cannot be used by sourcing its files:

```r
library(NADIA)          # correct
# source("R/Processing.R")   # does NOT work: the package is not a script bundle
```

Requires R >= 4.6. If NADIA is not installed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("NADIA")           # once released
# remotes::install_github("sciordia/NADIA")   # meanwhile
```

Many methods live in `Suggests:` and are only needed if selected. When one is
missing, NADIA names the package and how to install it — read that message
rather than switching method silently.

## The shortest analysis that is defensible

```r
library(NADIA)
data(nadia_dia)

res <- process_proteomics(nadia_dia, verbose = FALSE)

table(res$DEPs_results$Comparison, res$DEPs_results$Change)
```

`process_proteomics()` runs normalisation → optional batch correction →
imputation → differential abundance, and returns everything in memory. It
writes no files unless `export_dir` is given.

**But do not stop there, and do not present those defaults as a result.** The
defaults (`cycloess` + `combo`) are reasonable starting points, not answers.
Section *Choosing the methods* below is the part that matters.

## What you get back

`process_proteomics()` returns a list with six elements:

| Element | What it is |
|---|---|
| `se_proc` | `SummarizedExperiment`: one assay per stage (`raw`, `log2`, `<norm_method>`, `<imputed>`) |
| `DEPs_results` | One row per protein group **per comparison**; the statistics |
| `BoxPlot_Input` | Long format, one row per protein × sample × assay |
| `PCA_Input` | Long format, one row per protein × sample, with significance flags |
| `comparisons` | The contrasts run, as `"numerator-denominator"` |
| `parameters` | Every setting used, plus the recorded `call` |

`DEPs_results` columns: `Protein.IDs`, `Comparison`, `Gene.Names`, `logFC`,
`P.Value`, `adj.P.Val`, `Change`, `Assay`, and the four missingness columns
`MissGlobal`, `MissComp`, `MissCND1`, `MissCND2`.

## The workflow, and where to read more

Each stage has a reference file. Read the one you need; do not load them all.

1. **Import** — `reference/import.md`
   Four readers, one output contract. How the experimental design is declared,
   and what to check before going further.
2. **Choosing the methods** — `reference/method-choice.md`
   13 normalisation and 20 imputation methods. How to choose on evidence
   instead of by habit. **Read this one; it is where analyses go wrong.**
3. **Processing** — `reference/processing.md`
   `process_proteomics()` argument by argument, batch correction, limma vs
   limpa, and how to run the stages separately.
4. **Validation** — `reference/validation.md`
   `normalization_metrics()`, `imputation_metrics()`, and the spike-in
   benchmark. How to justify the pipeline you chose.
5. **Visualisation** — `reference/visualisation.md`
   Which figure answers which question, interactive and static.
6. **Export and archiving** — `reference/export.md`
   TSV, Parquet, the single-file `.nadia` archive, and interactive HTML.

`scripts/run_pipeline.R` is a parameterised end-to-end script to adapt rather
than to run blindly.

## Guardrails

These are the mistakes the package's own defaults and history invite. Check
them before reporting any result.

**Condition names must be syntactically valid R names.** No hyphens, no spaces,
no leading digits. `makeContrasts` builds an expression from them, so `"cond-1"`
or `"7d"` fails, sometimes confusingly. Rename at import.

**A comparison string is `"numerator-denominator"`.** `"B-A"` is B relative to
A. A positive `logFC` means higher in the numerator.

**`max_na_prop` is `NULL` by default, meaning no filtering.** It used to apply
to some imputation methods and not others. If you want proteins dropped by
missingness, pass a value explicitly; it then applies to every method alike.

**`halfmin` is `min - 1`, not `min / 2`.** The assay is log2, so halving an
intensity is subtracting 1. Do not "correct" it.

**`pattern_profiler_analysis()` clusters the assay named in
`DEPs_results$Assay`.** With the default `assay_name = NULL` it is resolved
from there, because that column records the assay the differential abundance
was computed on, and `assay_name` also filters those results — the two have to
agree. Pass it explicitly to cluster something else, or when the results cover
more than one assay, in which case NADIA refuses to guess. The names are in
`SummarizedExperiment::assayNames(res$se_proc)`.

**The four missingness columns are not interchangeable.** `MissGlobal` spans
every sample in the experiment; `MissComp` only the two conditions compared;
`MissCND1`/`MissCND2` each condition alone. A protein with `MissComp = 0` and
`MissGlobal = 60` is fully observed in this contrast — the missing values are
elsewhere. Do not filter on the wrong one.

**Batch correction with ComBat needs the biological variable.** Pass condition
in `batch_covariates`, or ComBat removes all batch-associated variance and can
erase the signal when condition and batch are confounded. NADIA warns when the
argument is `NULL`; do not ignore the warning.

**Pattern Profiler needs `Biobase` and `e1071` attached**, not merely
installed — Mfuzz calls them unqualified. `pattern_profiler_analysis()` handles
this itself, but a hand-rolled Mfuzz call will not.

**Never keep an open `.nadia` file inside iCloud Drive, Dropbox or OneDrive.**
Write and use it locally, close it, then copy the finished file across.

**A result is not a finding.** Imputation invents the values it fills. Before
reporting a protein as changed, look at its missingness columns: a large
`logFC` resting on imputed values in one condition is a presence/absence
observation, not a measured fold change. Say so.

## Reporting

When you present results, state the normalisation method, the imputation
method (and its MAR/MNAR stages if hybrid), the DE method, `alpha` and
`logFC_threshold`. Without them the numbers cannot be interpreted or
reproduced. `res$parameters` holds all of it, including the original call.
