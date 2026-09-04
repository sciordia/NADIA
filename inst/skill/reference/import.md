# Import: from a quantification report to `proteomics_data`

Four readers, one output contract. Everything downstream consumes any of them
unchanged, so this is the only stage that depends on the search engine.

| Report | Function | Shape |
|---|---|---|
| Spectronaut | `preprocess_spectronaut()` | **long**: one row per protein group per run |
| DIA-NN `report.pg_matrix.tsv` | `preprocess_diann()` | wide |
| TMT from Proteome Discoverer | `preprocess_tmt()` | wide |
| Label-free (DDA/DIA) from Proteome Discoverer | `preprocess_lfq()` | wide |

All four return an S3 object of class
`c("<type>_data", "proteomics_data", "list")` with exactly three elements:

- `metadata` — one row per sample
- `protein_id` — protein annotations and per-sample identification metrics
- `protein_quant` — the quantification matrix plus annotation columns

## Signatures

```r
preprocess_spectronaut(file_path, condition_order, export_dir = NULL,
                       agg_coverage_run = "max", agg_coverage_global = "max",
                       agg_mw = "max", agg_cscore_runwise = "mean",
                       timestamp_suffix = TRUE, verbose = TRUE)

preprocess_diann(file_path, condition_order, annot_path = NULL, export_dir = NULL, ...)
preprocess_tmt  (file_path, condition_order, annot_path = NULL, export_dir = NULL, ...)
preprocess_lfq  (file_path, annot_path = NULL, condition_order = NULL, export_dir = NULL, ...)
```

`condition_order` sets the factor level order, and the first level is the
reference. It decides more than the sign: it decides **which contrasts exist
and what they are called**. With `c("A", "B", "D")` you get `B-A`, `D-A`,
`D-B`; with `c("D", "B", "A")` you get `B-D`, `A-D`, `A-B` instead. The same
pair now has the opposite name and the opposite sign — `A-B` correlates −1 with
`B-A` — so downstream code that selects `comparisons = "B-A"` silently finds
nothing. Give the order deliberately and keep it fixed across a project.

`preprocess_lfq()` is the only reader where the argument is optional. Omitted,
the conditions fall in alphabetical order, which is a decision made for you
rather than by you: pass it anyway.

## How the design is declared in the three wide readers

Two routes, and only two. NADIA never guesses the design from a raw file name.

**Route 1 — column suffixes.** Columns named `Abundance: <condition>_<replicate>`:

```
Abundance: A_1   Abundance: A_2   Abundance: B_1   ...
```

Nothing else to supply. The replicate is the trailing `_<digits>`.

**Route 2 — an annotation sheet**, through `annot_path`. Two columns:

- `Column` — the sample identifier **exactly as it appears in the file**: the
  text after `Abundance: `, or the raw header for an unrenamed DIA-NN matrix.
- `Condition` — the condition it belongs to.

Matching is by name. It matches, or it aborts — there is no positional
fallback, on purpose: Proteome Discoverer does not write its column families in
a stable order, so a positional vector would silently relabel runs after a
re-export. Everything else in the sheet is ignored; **batch structure does not
go here** (see below).

A sheet that merely restates the suffixes produces the same data: `metadata`,
`protein_id` and `protein_quant` all come back `identical()`. The objects
themselves do not compare equal, because each carries provenance attributes
that record what was actually called — see below.

## Provenance travels with the object

Every `preprocess_*()` attaches two attributes. They are not list elements, so
the three-element contract is unchanged, and they are what `write_nadia()`
stores:

```r
attr(prep, "nadia_call")     # the call, as issued
attr(prep, "nadia_source")   # one row per input file: role, path, size, mtime, md5
```

`nadia_source` records the annotation sheet too when one was used, each with an
md5. Report them rather than describing the input in prose: they identify the
exact file the numbers came from. It also means `identical()` between two
objects built by different calls is `FALSE` even when the data match — compare
the three elements, or strip the attributes, if that is what you meant to test.

## Batch structure does not belong to the design

There is no batch column in a `proteomics_data` object. Batch is supplied later,
to `process_proteomics()`, through `covariate_df` + `batch_column`. Keep the
mapping (sample → batch) in your own table, and build it with the key NADIA
joins on:

```r
batch_df <- data.frame(
  Column = prep$metadata$Coding,        # <- the key. NOT R.FileName.
  Batch  = ifelse(prep$metadata$R.Replicate <= 4, "Mix1", "Mix2"))
```

`covariate_df` must have a column literally called `Column`, and its values are
matched against `metadata$Coding` — `"A_1"`, not the `"Abundance: A_1"` of
`R.FileName`. Get it wrong and the merge yields NAs; NADIA then aborts saying
the frame does not cover all samples.

## Verify the import before going on

Always, and report what you find:

```r
prep <- preprocess_spectronaut(path, condition_order = c("A", "B", "D"))

# 1. Did the design come out as intended?
table(prep$metadata$R.Condition)

# 2. How much is missing, before any imputation? Report the COUNT as well as
#    the percentage: a rounded percentage hides a handful of missing cells, and
#    a handful is the difference between "nothing to impute" and "something".
q <- as.matrix(prep$protein_quant[, grep("PG.Quantity", names(prep$protein_quant))])
n_missing <- sum(is.na(q) | q == 0)
sprintf("%d of %d cells (%.3f %%)", n_missing, length(q), 100 * n_missing / length(q))

# 3. Size
nrow(prep$protein_quant)
```

The `PG.Quantity_*` naming is shared by all four readers, so this snippet works
whichever one produced the object.

Zeros count as missing: NADIA converts them during normalisation. A report that
looks 0 % missing because the exporter wrote zeros is not complete data.

If the missingness is far from what the experiment should give (say 40 % on a
DIA run of four replicates), stop and ask about the export settings before
analysing. No imputation method rescues a bad export.

### When there is (almost) nothing missing

Some experiments arrive complete — the TMT example shipped with the package has
one missing cell in 48,000. Say so, and draw the consequence rather than going
through the motions:

- The imputation method is **irrelevant**: with nothing to fill, every method
  returns the input. Do not run `imputation_metrics()`; its simulation masks
  observed values, so it would measure recovery of an artificial problem that
  this dataset does not have.
- Normalisation still matters, and so does batch correction. Spend the effort
  there.
- Report it plainly: "no imputation was required (1 missing value in 48,000)"
  is a stronger statement than any imputation ranking.

## Example data shipped with the package

```r
data(nadia_dia)                                  # ready-made proteomics_data
system.file("extdata", package = "NADIA")        # the trimmed reports
```

`nadia_dia`: 3 conditions (A, B, D) × 4 replicates, 2,000 protein groups,
8.1 % missing. The `extdata` reports cover all four formats and share the same
conditions; the Spectronaut, DIA-NN and LFQ ones are the same twelve
injections searched three ways, which makes them useful for showing how
processing choices move missingness. The TMT report is a different experiment
and the only one with a real batch structure (two mixes, replicates 1-4 and
5-8, bridged by four `IS` channels).

The 2,000 proteins are a **random** sample, not the best-covered ones. That is
deliberate: selecting the best-quantified proteins produced a dataset with no
missing values at all, in which every imputation method is the identity.
