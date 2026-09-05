# Preprocess DIA-NN protein-group matrices

Converts a DIA-NN protein-group matrix (`report.pg_matrix.tsv`, wide TSV
format) into the same output structure as
[`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md):
three data.frames (`metadata`, `protein_id`, `protein_quant`) ready for
the downstream pipeline
([`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
and the associated modules).

DIA-NN names its intensity columns after the raw file, which says
nothing about the experimental design, so the design has to be declared.
There are two ways: rename the columns to
`Abundance: <Condition>_<Replicate>` as Proteome Discoverer writes them,
or pass a sample sheet through `annot_path`. The design is never
inferred from the raw-file name – a guess that happened to be right for
one export would be wrong and silent for the next.

A protein-group matrix carries only intensities. The identification
metrics the wide contract expects are filled from what DIA-NN does
report (`N.Sequences`, `N.Proteotypic.Sequences`), left as `NA` where it
reports nothing (the precursor counts), or omitted where nothing
downstream needs them (coverage and the Spectronaut C-score).

## Usage

``` r
preprocess_diann(
  file_path,
  condition_order,
  annot_path = NULL,
  export_dir = NULL,
  timestamp_suffix = TRUE,
  verbose = TRUE
)
```

## Arguments

- file_path:

  Path to the `report.pg_matrix.tsv` exported by DIA-NN.

- condition_order:

  Character vector with the order of the experimental conditions (e.g.
  `c("A", "B", "D")`). Only the runs whose condition is in this vector
  are kept, so it doubles as a filter. Conditions listed but absent from
  the report are dropped, with a warning, so that no empty factor level
  reaches the metadata.

- annot_path:

  Path to a sample sheet with columns `Column` and `Condition`, used to
  read the report exactly as DIA-NN wrote it: `Column` holds the run
  headers, raw-file names and all. If `NULL` (default), the columns are
  expected to follow the `Abundance:` convention. Any other column in
  the sheet is ignored – batch structure belongs in `covariate_df` of
  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md),
  see
  [`vignette("batch-correction")`](https://sciordia.github.io/NADIA/articles/batch-correction.md).

  The sheet is a file rather than a vector of labels because a vector
  would be positional, and nothing guarantees the column order of an
  export: a re-export that reordered them would relabel every run
  without any way to notice. Matching by name either matches or fails.

- export_dir:

  Directory to export the TSV files to. If `NULL` (default), no files
  are exported.

- timestamp_suffix:

  Logical. If `TRUE` (default), appends a timestamp to the names of the
  exported files.

- verbose:

  Logical. If `TRUE` (default), shows progress messages.

## Value

A list with class `c("diann_data", "proteomics_data", "list")`
containing:

- metadata:

  Data frame with one record per run, including the number of protein
  groups quantified in each

- protein_id:

  Data frame with the identification metrics per protein (wide format,
  the global metrics replicated per run)

- protein_quant:

  Data frame with the quantification metrics per protein (includes
  `PG.Quantity_<Coding>`)

The call that produced the object and a fingerprint of the file it was
read from (path, size, modification time and MD5) travel with it as the
attributes `nadia_call` and `nadia_source`. They are attributes rather
than list elements so that the three-element structure above is
unchanged;
[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)
records them as the provenance of an analysis.

## Examples

``` r
# A trimmed DIA-NN protein-group matrix ships with the package. Its run
# columns have been renamed to the "Abundance: <Condition>_<Replicate>"
# convention, so no further argument is needed.
report <- system.file("extdata", "nadia_diann_report.tsv.gz",
                      package = "NADIA")

diann <- preprocess_diann(report, condition_order = c("A", "B", "D"),
                          verbose = FALSE)
diann
#> Preprocessed DIA-NN data
#> ------------------------
#> Runs (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
table(diann$metadata$R.Condition)
#> 
#> A B D 
#> 4 4 4 

# Straight from DIA-NN, where the columns are still raw-file paths, the
# design is given instead by a sheet whose 'Column' holds those paths:
# preprocess_diann("report.pg_matrix.tsv",
#                  condition_order = c("A", "B", "D"),
#                  annot_path = "design.tsv")
```
