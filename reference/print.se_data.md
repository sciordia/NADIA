# Print a Summary of Data Imported from a SummarizedExperiment

Print a Summary of Data Imported from a SummarizedExperiment

## Usage

``` r
# S3 method for class 'se_data'
print(x, ...)
```

## Arguments

- x:

  An `se_data` object, as returned by
  [`preprocess_se()`](https://sciordia.github.io/NADIA/reference/preprocess_se.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the summary it prints.

## Examples

``` r
data(nadia_dia)
quant <- as.matrix(nadia_dia$protein_quant[, grep("^PG\\.Quantity_",
                                                  names(nadia_dia$protein_quant))])
rownames(quant) <- nadia_dia$protein_quant$PG.ProteinGroups
se <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(intensity = quant),
  colData = data.frame(Condition = nadia_dia$metadata$R.Condition))
print(preprocess_se(se, verbose = FALSE))
#> Protein intensities imported from a SummarizedExperiment
#> --------------------------------------------------------
#> Samples (metadata): 12 
#> Proteins: 2000 
#> 
#> Conditions: A, B, D 
```
