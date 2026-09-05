# Import confusion matrices from multiple benchmark result folders

Reads `benchmark_confusion_overall.tsv` files from subdirectories,
adding an `Assay` column derived from the folder name.

## Usage

``` r
import_confusion_results(
  results_dir,
  pattern = "benchmark_confusion_overall\\.tsv$",
  method_names = NULL,
  recursive = TRUE
)
```

## Arguments

- results_dir:

  Parent directory containing benchmark subfolders

- pattern:

  Regex pattern for the confusion filename

- method_names:

  Optional character vector of Assay names

- recursive:

  Search subdirectories recursively? (default: TRUE)

## Value

data.frame with columns: Assay, Comparison, N, TP, FP, TN, FN, plus
percentage columns

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
confusion_all <- import_confusion_results(root)
#> import_confusion_results: loaded 2 method(s), 2 comparison(s), 4 rows.
confusion_all[, c("Assay", "Comparison", "TP", "FP", "TN", "FN")]
#>      Assay Comparison  TP  FP   TN  FN
#> 1 cycloess        B-A 469  27 1227 274
#> 2 cycloess        D-A 693  60 1194  50
#> 3 quantile        B-A 494  35 1219 249
#> 4 quantile        D-A 709 917  337  34
```
