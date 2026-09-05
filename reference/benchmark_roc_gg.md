# ROC Curves by Comparison (ggplot2)

Generates ROC curves for each comparison using pROC::ggroc(). Requires
the pROC package.

## Usage

``` r
benchmark_roc_gg(
  classified_df,
  p_col = "adj.P.Val",
  comparisons = NULL,
  title = "ROC Curves by Comparison",
  zoom = FALSE,
  palette = "Set1"
)
```

## Arguments

- classified_df:

  Classified data frame (output of .classify_all_comparisons)

- p_col:

  P-value column name used as predictor score

- comparisons:

  Comparisons to include (NULL = all)

- title:

  Plot title

- zoom:

  If TRUE, zoom into the low-FPR region (x = 0 to 0.1)

- palette:

  RColorBrewer palette name (default: "Set1")

## Value

ggplot2 object or NULL if pROC is not available

## Examples

``` r
data(nadia_dia)
de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
if (requireNamespace("pROC", quietly = TRUE)) {
  bench <- benchmarking_proteomics(de, ev, verbose = FALSE)
  benchmark_roc_gg(bench$classified_df, comparisons = "B-A")
}
```
