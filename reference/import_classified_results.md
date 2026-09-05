# Import classified results from multiple benchmark result folders

Reads `benchmark_classified.tsv` files from subdirectories, adding an
`Assay` column derived from the folder name.

## Usage

``` r
import_classified_results(
  results_dir,
  pattern = "benchmark_classified\\.tsv$",
  method_names = NULL,
  recursive = TRUE
)
```

## Arguments

- results_dir:

  Parent directory containing benchmark subfolders

- pattern:

  Regex pattern for the classified filename

- method_names:

  Optional character vector of Assay names

- recursive:

  Search subdirectories recursively? (default: TRUE)

## Value

data.frame with columns from classified_df plus Assay

## Examples

``` r
data(nadia_dia)
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
root <- file.path(tempdir(), "nadia_bench")
for (nm in c("cycloess", "quantile")) {
  de <- process_proteomics(nadia_dia, norm_method = nm,
                           verbose = FALSE)$DEPs_results
  de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
  benchmarking_proteomics(de, ev, output_dir = file.path(root, nm),
                          verbose = FALSE)
}
classified_all <- import_classified_results(root)
#> import_classified_results: loaded 2 method(s), 2 comparison(s), 7988 rows.
table(classified_all$Assay, classified_all$classification)
#>           
#>              FN   FP   TN   TP
#>   cycloess  324   87 2421 1162
#>   quantile  283  952 1556 1203
```
