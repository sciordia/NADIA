# Compute dispersion metrics by Comparison x Species

Calculates MED, SD, CV, MAD, RCV, IQR, trimmed SD/CV for significant
proteins with correct direction.

## Usage

``` r
compute_dispersion_metrics(
  de_res,
  ev,
  alpha = 0.05,
  lfc_thr = 0,
  p_col = "adj.P.Val",
  comparisons = NULL,
  assay = NULL,
  trim = 0.1,
  species_df = NULL
)
```

## Arguments

- de_res:

  Data frame with DE results

- ev:

  Data frame with expected values

- alpha:

  Significance threshold (default: 0.05)

- lfc_thr:

  Log fold-change threshold (default: 0)

- p_col:

  P-value column name (default: "adj.P.Val")

- comparisons:

  Comparisons to include (NULL = all)

- assay:

  Assay to filter (NULL = all)

- trim:

  Trimming proportion for trimmed SD/CV (default: 0.1)

- species_df:

  Data frame with Protein.IDs and Species columns (optional). If NULL,
  de_res must already contain a Species column.

## Value

Data frame with dispersion metrics

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
disp <- compute_dispersion_metrics(de, ev)
disp[, c("Comparison", "Species", "expected_logFC", "MED", "MAD", "RCV")]
#>   Comparison Species expected_logFC     MED    MAD    RCV
#> 1        B-A   HUMAN             NA  0.1613 0.0793  49.13
#> 2        B-A   YEAST          -0.58 -0.5529 0.0582  10.53
#> 3        B-A   ECOLI           1.00  1.1209 0.1002   8.94
#> 4        D-A   HUMAN             NA  0.1301 0.2073 159.30
#> 5        D-A   YEAST          -3.30 -6.8168 1.1267  16.53
#> 6        D-A   ECOLI           2.00  2.1061 0.1326   6.30
```
