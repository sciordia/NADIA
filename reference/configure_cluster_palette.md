# Configure the colour palette for the clusters

Configure the colour palette for the clusters

## Usage

``` r
configure_cluster_palette(n_clusters, palette = NULL)
```

## Arguments

- n_clusters:

  Number of clusters

- palette:

  Colour palette: `NULL` for the default; `"pkg::name"` for a palette
  from paletteer; `"brewer:Name"`, or the bare name, for an RColorBrewer
  palette; or a vector of colours. A palette with fewer colours than
  clusters is interpolated.

## Value

Vector of colours named by cluster
