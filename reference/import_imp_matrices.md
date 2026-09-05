# Import imputed matrices into a SummarizedExperiment

Reads all TSV files matching a pattern from a directory, one assay per
file. Each TSV must have column 1 as protein identifier and remaining
columns as sample names matching the `sample_col` in the metadata file.

## Usage

``` r
import_imp_matrices(
  tsv_dir,
  metadata_path,
  pattern = "matrix_log2_.*\\.tsv$",
  method_names = NULL,
  condition_col = "Condition",
  sample_col = "Column"
)
```

## Arguments

- tsv_dir:

  Character. Directory containing the TSV files.

- metadata_path:

  Character. Path to a TSV or CSV file with at minimum columns
  `sample_col` and `condition_col`.

- pattern:

  Character. Regex to filter files in `tsv_dir`. Default:
  `"matrix_log2_.*\\.tsv$"`.

- method_names:

  Character vector. Assay names to assign (in order). NULL (default) =
  extracted from file names between `matrix_log2_` and `.tsv`.

- condition_col:

  Character. Column in metadata with condition/group labels.

- sample_col:

  Character. Column in metadata with sample names.

## Value

A SummarizedExperiment with one assay per TSV, colData from metadata,
and rowData containing `Protein.IDs`.

## Examples

``` r
# Round trip: write the imputed matrices out, then read them back in
data(nadia_dia)
res <- process_proteomics(nadia_dia, export_dir = tempdir(),
                          export_format = "tsv", verbose = FALSE)

meta <- file.path(tempdir(), "metadata.tsv")
utils::write.table(
  data.frame(Column = res$se_proc$Column, Condition = res$se_proc$Condition,
             Replicate = res$se_proc$Replicate),
  meta, sep = "\t", row.names = FALSE, quote = FALSE)

se <- import_imp_matrices(tsv_dir = tempdir(), metadata_path = meta)
#> import_imp_matrices: loaded 2 assay(s) -- cycloess, cycloess_Impseqrob_min
#>   Proteins: 1997 | Samples: 12
SummarizedExperiment::assayNames(se)
#> [1] "cycloess"               "cycloess_Impseqrob_min"
```
