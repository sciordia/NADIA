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

`condition_order` sets the factor level order, which decides the direction of
every contrast. Give it deliberately: `c("A", "B", "D")` makes A the reference.

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

A sheet that merely restates the suffixes gives a result `identical()` to using
no sheet.

## Batch structure does not belong to the design

There is no batch column in a `proteomics_data` object. Batch is supplied later,
to `process_proteomics()`, through `covariate_df` + `batch_column`. Keep the
mapping (sample → batch) in your own table.

## Verify the import before going on

Always, and report what you find:

```r
prep <- preprocess_spectronaut(path, condition_order = c("A", "B", "D"))

# 1. Did the design come out as intended?
table(prep$metadata$R.Condition)

# 2. How much is missing, before any imputation?
q <- prep$protein_quant[, grep("PG.Quantity", names(prep$protein_quant))]
round(100 * mean(is.na(as.matrix(q)) | as.matrix(q) == 0), 1)

# 3. Size
nrow(prep$protein_quant)
```

Zeros count as missing: NADIA converts them during normalisation. A report that
looks 0 % missing because the exporter wrote zeros is not complete data.

If the missingness is far from what the experiment should give (say 40 % on a
DIA run of four replicates), stop and ask about the export settings before
analysing. No imputation method rescues a bad export.

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
