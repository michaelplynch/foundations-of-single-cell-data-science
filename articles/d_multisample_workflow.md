# scRNA-seq multi-sample analysis

``` r


library(dittoSeq)
library(scuttle)
library(muscData)
```

## Multi-sample Single-cell RNA-seq in Bioconductor

Multi-sample scRNA-seq analysis broadly follows the same analysis step
as for single-sample, but there are some key differences. The process is
well described in
[OSCA](https://bioconductor.org/books/release/OSCA.multisample/).

Below, a demo multi-sample dataset is loaded, containing single-cell
data from 8 matched lupus patients before and after IFN stimulation,
which can be used for self-directed learning.

### Data loading

Data used are from Kang HM, et al. (2018). Multiplexed droplet
single-cell RNA-sequencing using natural genetic variation. *Nature
Biotechnology* 36:89–94.
[doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042). PBMCs from 8
lupus patients, control and IFN-β stimulated, distributed via the
muscData Bioconductor package.

``` r


# find out more about the dataset
?Kang18_8vs8

sce<-Kang18_8vs8()
#> see ?muscData and browseVignettes('muscData') for documentation
#> loading from cache
```

### Analysis

The rest is up to you :)
