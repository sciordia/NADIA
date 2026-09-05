# Import normalized matrices into a SummarizedExperiment

Reads all TSV files matching a pattern from a directory, one assay per
file. Each TSV must have column 1 named `ProteinGroups` and remaining
columns as sample names matching the `sample_col` in the metadata file.

## Usage

``` r
import_norm_matrices(
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

  Character. Directory containing the normalized TSV files.

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
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se <- nm_run_normalizations(se, methods = c("cycloess", "MAD"), verbose = FALSE)

# Write the normalized assays as TSVs, then read them back into a new SE
tsv_dir <- file.path(tempdir(), "nm_matrices")
dir.create(tsv_dir, showWarnings = FALSE)
for (a in c("cycloess", "MAD")) {
  m <- SummarizedExperiment::assay(se, a)
  utils::write.table(
    data.frame(ProteinGroups = rownames(m), m, check.names = FALSE),
    file.path(tsv_dir, paste0("matrix_log2_", a, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE)
}
meta_path <- file.path(tsv_dir, "metadata.tsv")
utils::write.table(as.data.frame(SummarizedExperiment::colData(se)), meta_path,
                   sep = "\t", row.names = FALSE, quote = FALSE)

se_nm <- import_norm_matrices(tsv_dir, meta_path)
#> import_norm_matrices: loaded 2 assay(s) -- MAD, cycloess
#>   Proteins: 1997 | Samples: 12
SummarizedExperiment::assayNames(se_nm)
#> [1] "MAD"      "cycloess"
```
