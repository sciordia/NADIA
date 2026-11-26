# Rename rownames using IDs column

Rename rownames using IDs column

## Usage

``` r
.rename_rownames_from_rd(x_df, rd, id_col = "IDs")
```

## Arguments

- x_df:

  Matrix or data frame with rownames as ProteinGroups

- rd:

  rowData data frame with IDs column

- id_col:

  Name of IDs column (default: "IDs")

## Value

Data frame with renamed rownames
