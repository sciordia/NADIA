# Build a PCA scores data frame from data in long format

Build a PCA scores data frame from data in long format

## Usage

``` r
build_pca_scores(
  pca_input,
  mode = c("all", "any", "specific"),
  alpha = 0.05,
  comparison = NULL,
  subset_label = NULL,
  center = TRUE,
  scale. = TRUE,
  filter_samples_to_comparison = FALSE,
  cond_col = "Condition"
)
```

## Arguments

- pca_input:

  Data frame in long format with columns:

  - SampleID: Sample identifier

  - FeatureID: Protein/feature identifier

  - Intensity: Intensity value (log2)

  - Condition: Experimental condition

  - Replicate: Replicate number (optional)

  - sig_any: Logical flagging significance in any comparison (for
    mode="any")

  - adjP\_\*: Adjusted p-value columns, one per comparison (for
    mode="specific")

- mode:

  Protein filtering mode: "all", "any", or "specific"

- alpha:

  Significance threshold for mode "specific" (default: 0.05)

- comparison:

  Name of the comparison for mode "specific" (e.g. "B-A")

- subset_label:

  Custom label for the subset (optional)

- center:

  Center the data before PCA (default: TRUE)

- scale.:

  Scale the data before PCA (default: TRUE)

- filter_samples_to_comparison:

  Restrict samples to the conditions involved in the specific comparison
  (default: FALSE)

- cond_col:

  Name of the condition column (default: "Condition")

## Value

Data frame with columns: SampleID, PC1, PC2, PC1_Perc, PC2_Perc, Subset,
Condition, Replicate

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# Scores and the variance explained by the first two components
scores <- build_pca_scores(res$PCA_Input, mode = "all")
head(scores)
#>   SampleID        PC1        PC2 PC1_Perc PC2_Perc       Subset Condition
#> 1      A_1 -31.330282   3.479060    37.59    11.33 All proteins         A
#> 2      A_2 -32.180712  30.462243    37.59    11.33 All proteins         A
#> 3      A_3 -26.921997  -7.797913    37.59    11.33 All proteins         A
#> 4      A_4 -24.416436 -20.914206    37.59    11.33 All proteins         A
#> 5      B_1  -8.655200   9.415768    37.59    11.33 All proteins         B
#> 6      B_2  -6.694873   6.716063    37.59    11.33 All proteins         B
#>   Replicate
#> 1         1
#> 2         2
#> 3         3
#> 4         4
#> 5         1
#> 6         2
unique(scores[c("PC1_Perc", "PC2_Perc")])
#>   PC1_Perc PC2_Perc
#> 1    37.59    11.33
```
