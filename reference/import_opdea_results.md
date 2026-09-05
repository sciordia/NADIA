# Import OpDEA metrics from multiple benchmark result folders

Reads `benchmark_opdea_metrics.tsv` files from subdirectories, adding an
`Assay` column derived from the folder name.

## Usage

``` r
import_opdea_results(
  results_dir,
  pattern = "benchmark_opdea_metrics\\.tsv$",
  method_names = NULL,
  recursive = TRUE
)
```

## Arguments

- results_dir:

  Parent directory containing benchmark subfolders

- pattern:

  Regex pattern for the opdea metrics filename (default:
  `"benchmark_opdea_metrics\.tsv$"`)

- method_names:

  Optional character vector of Assay names. If NULL, names are derived
  from the parent folder of each file.

- recursive:

  Logical. Search subdirectories recursively? (default: TRUE)

## Value

data.frame in long format with columns: Assay, Comparison, nMCC, G_mean,
pAUC_001, pAUC_005, pAUC_010

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
opdea_all <- import_opdea_results(root)
#> import_opdea_results: loaded 2 method(s), 2 comparison(s), 4 rows.
head(opdea_all)
#>   Comparison   nMCC G_mean pAUC_001 pAUC_005 pAUC_010    Assay
#> 1        B-A 0.8411 0.7859   0.7677   0.8071   0.8211 cycloess
#> 2        D-A 0.9412 0.9424   0.9125   0.9505   0.9614 cycloess
#> 3        B-A 0.8488 0.8040   0.7173   0.8046   0.8305 quantile
#> 4        D-A 0.6386 0.5064   0.7150   0.7554   0.7695 quantile
```
