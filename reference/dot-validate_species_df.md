# Validate and de-duplicate the species mapping

`species_df` is joined to the DE results on `Protein.IDs`, so a repeated
identifier multiplies the rows of every comparison it appears in and
inflates the counts the metrics are built from. The mapping is therefore
reduced to one row per protein before the join: exact repeats are
dropped, and an identifier claimed by two different species is an error,
because there is no way to know which side of the truth table the
protein belongs to.

## Usage

``` r
.validate_species_df(species_df)
```

## Arguments

- species_df:

  Data frame with Protein.IDs and Species columns

## Value

The mapping reduced to the two columns, one row per protein
