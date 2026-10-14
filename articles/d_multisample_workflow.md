# scRNA-seq multi-sample analysis

``` r


library(dittoSeq)
library(scuttle)
library(muscData)
```

## Data loading

``` r


sce<-Kang18_8vs8()
#> see ?muscData and browseVignettes('muscData') for documentation
#> loading from cache

citation("muscData")
#> To cite package 'muscData' in publications use:
#> 
#>   Crowell H (2026). _muscData: Multi-sample multi-group scRNA-seq
#>   data_. doi:10.18129/B9.bioc.muscData
#>   <https://doi.org/10.18129/B9.bioc.muscData>. R package version
#>   1.26.0, <https://bioconductor.org/packages/muscData>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {muscData: Multi-sample multi-group scRNA-seq data},
#>     author = {Helena L. Crowell},
#>     year = {2026},
#>     note = {R package version 1.26.0},
#>     url = {https://bioconductor.org/packages/muscData},
#>     doi = {10.18129/B9.bioc.muscData},
#>   }
```

## Analysis

The rest is up to you :)
