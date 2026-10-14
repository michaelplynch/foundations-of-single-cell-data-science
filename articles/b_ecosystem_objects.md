# Bioconductor Ecosystem and the SingleCellExperiment Object

## Bioconductor

[Bioconductor](https://www.bioconductor.org/) is a project which aims to
develop and share open-source software for precise and repeatable
analysis of biological data. They foster an inclusive and collaborative
community of developers and data scientists.

## Why Bioconductor?

- Analysis workflows are long and complex, an ecosystem for R packages
  for genomics data analysis facilitates reuse.
- Common data structures allow interoperability between packages.
- Package build system automates testing, ensures compatibility, and
  improves reproducibility.
- Community events and forums improve collaboration.

You can:

- Use high-quality, interoperable, and well-documented Bioconductor R
  packages.
- Develop and
  [contribute](https://contributions.bioconductor.org/bioconductor-package-submissions.html)
  your own R packages for genomics analysis.

## Accessing Bioconductor packages

You might typically install a package from CRAN through
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html).
Bioconductor packages are hosted elsewhere and have a different
mechanism for installing.

First, we’ll install BiocManager from CRAN, and we can then use
BiocManager to install packages from Bioconductor.

``` r


## not required as packages installed through workshop package installation.
install.packages('BiocManager')
BiocManager::install('SingleCellExperiment')
```

You can then load packages as normal.

``` r


library(SingleCellExperiment)
library(scrapper)
library(TENxPBMCData)
```

As well as hosting, testing, and distributing packages, Bioconductor
also ships some core infrastructure which improves interoperability
between packages. These usually come in the form of data structures for
storing different types of data. In this workshop, we’ll use the
SingleCellExperiment class.

We begin by loading in a public dataset from the TENxPBMCData
Bioconductor package. The dataset contains approx. 4,000 cells from
peripheral blood mononuclear cells (PBMCs). These are mostly immune
cells.

Data used are from the *4k PBMCs from a Healthy Donor* dataset released
by 10x Genomics (2017), distributed via the TENxPBMCData Bioconductor
package.

``` r


sce<-TENxPBMCData::TENxPBMCData('pbmc4k')
#> see ?TENxPBMCData and browseVignettes('TENxPBMCData') for documentation
#> loading from cache

citation("TENxPBMCData")
#> To cite package 'TENxPBMCData' in publications use:
#> 
#>   Hansen KD, Risso D, Hicks S (2026). _TENxPBMCData: PBMC data from 10X
#>   Genomics_. doi:10.18129/B9.bioc.TENxPBMCData
#>   <https://doi.org/10.18129/B9.bioc.TENxPBMCData>. R package version
#>   1.30.0, <https://bioconductor.org/packages/TENxPBMCData>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {TENxPBMCData: PBMC data from 10X Genomics},
#>     author = {Kasper D. Hansen and Davide Risso and Stephanie Hicks},
#>     year = {2026},
#>     note = {R package version 1.30.0},
#>     url = {https://bioconductor.org/packages/TENxPBMCData},
#>     doi = {10.18129/B9.bioc.TENxPBMCData},
#>   }
```

The SingleCellExperiment object keeps the different types of data
relevant to our analysis neatly organised in one object. Each data point
is associated with a feature and an observation, but we might also want
to keep other information, for example:

- each feature (gene) can be named using different conventions, or have
  different aliases.
- each observation (cell) might need its patient ID, sequencing batch,
  quality metrics or cluster ID attached.

``` r


class(sce)
#> [1] "SingleCellExperiment"
#> attr(,"package")
#> [1] "SingleCellExperiment"

sce
#> class: SingleCellExperiment 
#> dim: 33694 4340 
#> metadata(0):
#> assays(1): counts
#> rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#>   ENSG00000268674
#> rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
#> colnames: NULL
#> colData names(11): Sample Barcode ... Individual Date_published
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

counts(sce)[1:4,1:4]
#> <4 x 4> sparse DelayedMatrix object of type "integer":
#>                 [,1] [,2] [,3] [,4]
#> ENSG00000243485    0    0    0    0
#> ENSG00000237613    0    0    0    0
#> ENSG00000186092    0    0    0    0
#> ENSG00000238009    0    0    0    0

head(colData(sce))
#> DataFrame with 6 rows and 11 columns
#>        Sample            Barcode         Sequence   Library Cell_ranger_version
#>   <character>        <character>      <character> <integer>         <character>
#> 1      pbmc4k AAACCTGAGAAGGCCT-1 AAACCTGAGAAGGCCT         1                v2.1
#> 2      pbmc4k AAACCTGAGACAGACC-1 AAACCTGAGACAGACC         1                v2.1
#> 3      pbmc4k AAACCTGAGATAGTCA-1 AAACCTGAGATAGTCA         1                v2.1
#> 4      pbmc4k AAACCTGAGCGCCTCA-1 AAACCTGAGCGCCTCA         1                v2.1
#> 5      pbmc4k AAACCTGAGGCATGGT-1 AAACCTGAGGCATGGT         1                v2.1
#> 6      pbmc4k AAACCTGCAAGGTTCT-1 AAACCTGCAAGGTTCT         1                v2.1
#>   Tissue_status Barcode_type   Chemistry Sequence_platform    Individual
#>     <character>  <character> <character>       <character>   <character>
#> 1            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#> 2            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#> 3            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#> 4            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#> 5            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#> 6            NA     Chromium Chromium_v2         Hiseq4000 HealthyDonor1
#>   Date_published
#>      <character>
#> 1     2017-11-08
#> 2     2017-11-08
#> 3     2017-11-08
#> 4     2017-11-08
#> 5     2017-11-08
#> 6     2017-11-08

head(rowData(sce))
#> DataFrame with 6 rows and 3 columns
#>                      ENSEMBL_ID   Symbol_TENx       Symbol
#>                     <character>   <character>  <character>
#> ENSG00000243485 ENSG00000243485  RP11-34P13.3           NA
#> ENSG00000237613 ENSG00000237613       FAM138A      FAM138A
#> ENSG00000186092 ENSG00000186092         OR4F5        OR4F5
#> ENSG00000238009 ENSG00000238009  RP11-34P13.7 LOC100996442
#> ENSG00000239945 ENSG00000239945  RP11-34P13.8           NA
#> ENSG00000239906 ENSG00000239906 RP11-34P13.14           NA
```

For convenience, column data (data describing the cells) can be accessed
using the ‘\$’ operator.

``` r


table(sce$Individual)
#> 
#> HealthyDonor1 
#>          4340

table(sce$Individual)
#> 
#> HealthyDonor1 
#>          4340
```

We can also store transformations of the raw data, e.g. normalised
counts.

``` r


assays(sce)
#> List of length 1
#> names(1): counts
sce<-scrapper::normalizeRnaCounts.se(sce)
#> using unknown matrix fallback for 'HDF5ArraySeed'
assays(sce)
#> List of length 2
#> names(2): counts logcounts
```

Aaaand we can also store embeddings of the data, such as from PCA, TSNE
or UMAP.

``` r


reducedDimNames(sce) ## empty for now but we'll add to it later.
#> character(0)
```

Finally, we can store alternative experiments, or modalities, which were
measured on the same cells. For example,
[CITE-seq](https://emea.illumina.com/techniques/sequencing/rna-sequencing/cite-seq.html)
data measures gene and surface protein expression simultaneously.

Why might we want to treat these separately rather than appending them
to the main counts?

``` r


altExpNames(sce) ## also expected to be empty for this workshop.
#> character(0)
```

## Exercise

- Can you extract basic summary statistics from the metadata?
- Can you extract the counts data?
- Can you describe the counts data (e.g. distribution, sparsity)?

``` r


sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Etc/UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] TENxPBMCData_1.30.0         HDF5Array_1.40.0           
#>  [3] h5mread_1.4.1               rhdf5_2.56.1               
#>  [5] DelayedArray_0.38.2         SparseArray_1.12.3         
#>  [7] S4Arrays_1.12.1             abind_1.4-8                
#>  [9] Matrix_1.7-6                scrapper_1.6.3             
#> [11] SingleCellExperiment_1.34.0 SummarizedExperiment_1.42.0
#> [13] Biobase_2.72.0              GenomicRanges_1.64.0       
#> [15] Seqinfo_1.2.0               IRanges_2.46.0             
#> [17] S4Vectors_0.50.3            BiocGenerics_0.58.1        
#> [19] generics_0.1.4              MatrixGenerics_1.24.0      
#> [21] matrixStats_1.5.0           BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1     dplyr_1.2.1          blob_1.3.0          
#>  [4] Biostrings_2.80.2    filelock_1.0.3       fastmap_1.2.0       
#>  [7] BiocFileCache_3.2.0  digest_0.6.39        lifecycle_1.0.5     
#> [10] KEGGREST_1.52.2      RSQLite_3.53.3       magrittr_2.0.5      
#> [13] compiler_4.6.1       rlang_1.3.0          sass_0.4.10         
#> [16] tools_4.6.1          yaml_2.3.12          knitr_1.52          
#> [19] htmlwidgets_1.6.4    bit_4.6.0            curl_8.0.0          
#> [22] withr_3.0.3          purrr_1.2.2          desc_1.4.3          
#> [25] grid_4.6.1           ExperimentHub_3.2.2  beachmat_2.28.0     
#> [28] Rhdf5lib_2.0.0       cli_3.6.6            crayon_1.5.3        
#> [31] rmarkdown_2.32       ragg_1.5.2           otel_0.2.0          
#> [34] httr_1.4.9           DBI_1.3.0            cachem_1.1.0        
#> [37] parallel_4.6.1       AnnotationDbi_1.74.0 BiocManager_1.30.27 
#> [40] XVector_0.52.0       vctrs_0.7.3          jsonlite_2.0.0      
#> [43] bookdown_0.48        BiocNeighbors_2.6.0  bit64_4.8.6         
#> [46] systemfonts_1.3.2    jquerylib_0.1.4      glue_1.8.1          
#> [49] pkgdown_2.2.1        BiocVersion_3.23.1   tibble_3.3.1        
#> [52] pillar_1.11.1        rappdirs_0.3.4       htmltools_0.5.9     
#> [55] rhdf5filters_1.24.1  R6_2.6.1             dbplyr_2.6.0        
#> [58] httr2_1.3.0          textshaping_1.0.5    evaluate_1.0.5      
#> [61] lattice_0.22-9       AnnotationHub_4.2.2  png_0.1-9           
#> [64] memoise_2.0.1        bslib_0.12.0         Rcpp_1.1.2          
#> [67] xfun_0.61            fs_2.1.0             pkgconfig_2.0.3
```
