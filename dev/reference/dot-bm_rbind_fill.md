# Row-bind data.frames with heterogeneous columns

Combines a list of data.frames by the union of their columns, filling
missing columns with NA. Unlike `do.call(rbind, ...)`, this tolerates
methods whose exported tables carry different column sets (e.g. an
external tool like Proteome Discoverer whose classified table keeps its
own de_res columns).

## Usage

``` r
.bm_rbind_fill(df_list)
```
