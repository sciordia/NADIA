# Preprocess Proteome Discoverer TMT exports

Converts a Proteome Discoverer protein export (wide TSV format, with
`Abundance: <Condition>_<Replicate>` columns) into the same output
structure as
[`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md):
three data.frames (`metadata`, `protein_id`, `protein_quant`) ready for
the downstream pipeline
([`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
and the associated modules).

The global PD identification metrics (`# PSMs`, `# Peptides`,
`# Unique Peptides`, `Coverage [%]`) are replicated into per-channel
columns to fit the wide Spectronaut contract.

## Usage

``` r
preprocess_tmt(
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

  Path to the TSV exported from Proteome Discoverer.

- condition_order:

  Character vector with the order of the experimental conditions (e.g.
  `c("A","B","D","IS")`). Only the channels whose condition is in this
  vector are kept – handy for excluding Internal Standards by omitting
  `"IS"`. Conditions listed but absent from the report are dropped, with
  a warning, so that no empty factor level reaches the metadata.

- annot_path:

  Path to a sample sheet with columns `Column` and `Condition`, for
  reports whose channels are not named `<Condition>_<Replicate>` – an
  export left with the TMT tags, `126`, `127N` and so on. `Column` holds
  the text after `Abundance: `. If `NULL` (default), the design comes
  from the suffixes. Any other column in the sheet is ignored – batch
  structure belongs in `covariate_df` of
  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md),
  see
  [`vignette("batch-correction")`](https://sciordia.github.io/NADIA/articles/batch-correction.md).

- export_dir:

  Directory to export the TSV files to. If `NULL` (default), no files
  are exported.

- timestamp_suffix:

  Logical. If `TRUE` (default), appends a timestamp to the names of the
  exported files.

- verbose:

  Logical. If `TRUE` (default), shows progress messages.

## Value

A list with class `c("tmt_data", "proteomics_data", "list")` containing:

- metadata:

  Data frame with one record per channel/sample

- protein_id:

  Data frame with the identification metrics per protein (wide format,
  global metrics replicated per channel)

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
# A trimmed Proteome Discoverer TMTpro report ships with the package: three
# conditions of eight replicates, split across two TMT mixes (replicates 1-4
# are the first mix, 5-8 the second) and bridged by four "IS" channels.
report <- system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA")

# condition_order is required, and doubles as a filter: the report also holds
# an "IS" channel (internal standards) that is left out by not listing it
tmt <- preprocess_tmt(report,
                      condition_order = c("A", "B", "D"),
                      verbose = FALSE)
tmt
#> Preprocessed TMT (Proteome Discoverer) data
#> ------------------------------------------
#> Channels (metadata): 24 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
table(tmt$metadata$R.Condition)
#> 
#> A B D 
#> 8 8 8 
```
