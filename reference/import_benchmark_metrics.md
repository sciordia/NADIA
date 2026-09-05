# Import benchmark metrics from multiple benchmark result folders

Reads `benchmark_metrics.tsv` files from subdirectories, adding an
`Assay` column derived from the folder name.

## Usage

``` r
import_benchmark_metrics(
  results_dir,
  pattern = "benchmark_metrics\\.tsv$",
  method_names = NULL,
  recursive = TRUE
)
```

## Arguments

- results_dir:

  Parent directory containing benchmark subfolders

- pattern:

  Regex pattern for the metrics filename

- method_names:

  Optional character vector of Assay names

- recursive:

  Search subdirectories recursively? (default: TRUE)

## Value

data.frame with columns: Assay, Comparison, TP, FP, TN, FN, Sensitivity,
Specificity, Precision, NPV, Accuracy, F1, MCC, AUC, Performance

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
bench_all <- import_benchmark_metrics(root)
#> import_benchmark_metrics: loaded 2 method(s), 2 comparison(s), 4 rows.
bench_all[, c("Assay", "Comparison", "Sensitivity", "Specificity", "F1")]
#>      Assay Comparison Sensitivity Specificity     F1
#> 1 cycloess        B-A      0.6312      0.9785 0.7571
#> 2 cycloess        D-A      0.9327      0.9522 0.9265
#> 3 quantile        B-A      0.6649      0.9721 0.7767
#> 4 quantile        D-A      0.9542      0.2687 0.5986
```
