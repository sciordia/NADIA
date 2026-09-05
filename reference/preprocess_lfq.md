# Preprocess Proteome Discoverer LFQ exports

Converts a Proteome Discoverer LFQ protein export (wide TSV format) into
the same output structure as
[`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md)
/
[`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md):
three data.frames (`metadata`, `protein_id`, `protein_quant`) ready for
the downstream pipeline
([`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)).

The experimental design comes from the
`Abundance: <Condition>_<Replicate>` column suffixes, exactly as in
[`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md)
and
[`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md).
For a report whose columns are not named that way, `annot_path` takes a
sample sheet instead.

The intensities (`Abundance: <sample>`) and the three per-sample metric
families (`Score Mascot`, `# PSMs (by Search Engine)`,
`# Peptides (by Search Engine)`) are mapped BY SAMPLE NAME, never by
position: Proteome Discoverer writes the four families in different
orders.

## Usage

``` r
preprocess_lfq(
  file_path,
  annot_path = NULL,
  condition_order = NULL,
  export_dir = NULL,
  timestamp_suffix = TRUE,
  verbose = TRUE
)
```

## Arguments

- file_path:

  Path to the data TSV exported from Proteome Discoverer.

- annot_path:

  Path to a sample sheet with columns `Column` and `Condition`, where
  `Column` holds the text after `Abundance: `. If `NULL` (default), the
  design comes from the suffixes. Any other column in the sheet is
  ignored – batch structure belongs in `covariate_df` of
  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md),
  see
  [`vignette("batch-correction")`](https://sciordia.github.io/NADIA/articles/batch-correction.md).

- condition_order:

  Character vector with the order of the conditions. If `NULL`
  (default), it is derived from the data (order of appearance). If
  supplied, it fixes the factor levels and discards samples whose
  condition is not in the list; conditions listed but absent are dropped
  too, with a warning, so that no empty factor level reaches the
  metadata.

- export_dir:

  Directory to export the TSVs to. `NULL` (default) = no export.

- timestamp_suffix:

  Logical. If `TRUE` (default), appends a timestamp to the exported
  files.

- verbose:

  Logical. If `TRUE` (default), shows progress messages.

## Value

A list with class `c("lfq_data", "proteomics_data", "list")` containing
`metadata`, `protein_id` and `protein_quant`.

The call that produced the object and a fingerprint of the files it was
read from – the report, and the sample sheet when there is one – (path,
size, modification time and MD5) travel with it as the attributes
`nadia_call` and `nadia_source`. They are attributes rather than list
elements so that the three-element structure above is unchanged;
[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)
records them as the provenance of an analysis.

## Examples

``` r
# A trimmed Proteome Discoverer LFQ report ships with the package. Its columns
# follow the "Abundance: <Condition>_<Replicate>" convention, so the design
# needs no further argument.
report <- system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA")

lfq <- preprocess_lfq(report, verbose = FALSE)
lfq
#> Preprocessed LFQ (Proteome Discoverer) data
#> ------------------------------------------
#> Samples (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
lfq$metadata[, c("Coding", "R.Condition", "R.Replicate")]
#>     Coding R.Condition R.Replicate
#> A_1    A_1           A           1
#> A_2    A_2           A           2
#> A_3    A_3           A           3
#> A_4    A_4           A           4
#> B_1    B_1           B           1
#> B_2    B_2           B           2
#> B_3    B_3           B           3
#> B_4    B_4           B           4
#> D_1    D_1           D           1
#> D_2    D_2           D           2
#> D_3    D_3           D           3
#> D_4    D_4           D           4

# The same design, declared in a sheet instead. Use this route when the
# columns are not named by the convention; here it is the same experiment, so
# it produces the same sample table.
annot <- system.file("extdata", "nadia_lfq_annotation.tsv", package = "NADIA")
from_sheet <- preprocess_lfq(report, annot_path = annot, verbose = FALSE)
identical(from_sheet$metadata, lfq$metadata)
#> [1] TRUE
```
