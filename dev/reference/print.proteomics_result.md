# Print method for proteomics_result

Print method for proteomics_result

## Usage

``` r
# S3 method for class 'proteomics_result'
print(x, ...)
```

## Arguments

- x:

  proteomics_result object

- ...:

  Additional arguments (ignored)

## Value

`x`, invisibly. Called for the summary it prints: the assays of the
SummarizedExperiment, the comparisons performed, the number of
differential proteins and the parameters used.

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)
print(res)
#> === Proteomics Processing Result ===
#> 
#> SummarizedExperiment:
#>   - Proteins: 1997 
#>   - Samples: 12 
#>   - Assays: raw, log2, cycloess, Impseqrob_min 
#>   - Conditions: A, B, D 
#> 
#> Differential Results:
#>   - Total rows: 5991 
#>   - Comparisons: B-A, D-A, D-B 
#>     B-A: Up=325, Down=172
#>     D-A: Up=374, Down=383
#>     D-B: Up=286, Down=382
#> 
#> Parameters:
#>   - Normalization: cycloess 
#>     - Cyclic Loess method: fast 
#>     - Cyclic Loess iterations: 3 
#>     - Cyclic Loess span: 0.7 
#>   - Imputation: combo 
#>     - MAR method: Impseqrob 
#>     - MNAR method: min 
#>   - Alpha: 0.05 
#>   - logFC threshold: 0 
#>   - Output directory: (no export) 
```
