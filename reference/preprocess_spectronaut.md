# Preprocess Spectronaut reports

Converts a Spectronaut report (long TSV format) into three structured
tables: run metadata, protein identification and protein quantification.

The report has to carry the right columns, and rather than leave that to
be reconstructed by hand the package ships the Spectronaut export schema
that produces them:

    system.file("extdata", "NADIA_Report.rs", package = "NADIA")

Import that file into Spectronaut as a report schema and export with it.
The result has the shape of `nadia_dia_report.tsv.gz`, the example that
ships alongside it.

## Usage

``` r
preprocess_spectronaut(
  file_path,
  condition_order,
  export_dir = NULL,
  agg_coverage_run = c("max", "mean", "median", "min"),
  agg_coverage_global = c("max", "mean", "median", "min"),
  agg_mw = c("max", "mean", "median", "min"),
  agg_cscore_runwise = c("mean", "max", "median", "min"),
  timestamp_suffix = TRUE,
  verbose = TRUE
)
```

## Arguments

- file_path:

  Path to the Spectronaut TSV file.

- condition_order:

  Character vector that selects the experimental conditions to retain
  and sets their order in subsequent analysis stages (e.g.
  `c("Control", "Treated")`). Runs whose condition is not listed are
  dropped; conditions listed but absent from the report are dropped too,
  with a warning, so that no empty factor level reaches the metadata.

- export_dir:

  Directory to export the TSV files to. If `NULL` (default), no files
  are exported.

- agg_coverage_run:

  Aggregation method for the per-sample PG.Coverage. Options: "max"
  (default), "mean", "median", "min".

- agg_coverage_global:

  Aggregation method for PG.Coverage.Global. Options: "max" (default),
  "mean", "median", "min".

- agg_mw:

  Aggregation method for PG.MolecularWeight. Options: "max" (default),
  "mean", "median", "min".

- agg_cscore_runwise:

  Aggregation method for PG.Cscore.RunWise. Options: "mean" (default),
  "max", "median", "min".

- timestamp_suffix:

  Logical. If `TRUE` (default), appends a timestamp to the names of the
  exported files.

- verbose:

  Logical. If `TRUE` (default), shows progress messages.

## Value

A list with class `spectronaut_data` containing:

- metadata:

  Data frame with the run information (1 row per sample)

- protein_id:

  Data frame with the identification metrics per protein

- protein_quant:

  Data frame with the quantification metrics per protein

The call that produced the object and a fingerprint of the file it was
read from (path, size, modification time and MD5) travel with it as the
attributes `nadia_call` and `nadia_source`. They are attributes rather
than list elements so that the three-element structure above is
unchanged;
[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)
records them as the provenance of an analysis.

## Examples

``` r
# A trimmed Spectronaut report ships with the package
report <- system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA")

prep <- preprocess_spectronaut(report,
                               condition_order = c("A", "B", "D"),
                               verbose = FALSE)
prep
#> Preprocessed Spectronaut data
#> -----------------------------
#> Runs (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
dim(prep$protein_quant)
#> [1] 2000   44

# Nothing is written to disk unless export_dir is given
prep <- preprocess_spectronaut(report,
                               condition_order = c("A", "B", "D"),
                               export_dir = tempdir(),
                               verbose = FALSE)
```
