# The SingleCellExperiment Object and EDA

``` r


library(SingleCellExperiment)
library(scuttle)
library(muscData)
```

## Bioconductor

[Bioconductor](https://www.bioconductor.org/) is a project which aims to
develop and share open-source software for precise and repeatable
analysis of biological data. They foster an inclusive and collaborative
community of developers and data scientists.

## Why Bioconductor?

- Analysis workflows are long and complex, an ecosystem for R packages
  for genomics data analysis facilitates reuse.
- Common data structures allow interoperability between packages.
- Package build system automates testing, ensures compatability, and
  improves reprehensibility.
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


install.packages('BiocManager')
BiocManager::install('SingleCellExperiment')
```

As well as hosting, testing, and distributing packages, Bioconductor
also ships some core infrastructure which improves interoperability
between packages. These usually come in the form of data structures for
storing different types of data. In this workshop, we’ll use the
SingleCellExperiment class.

``` r


sce<-Kang18_8vs8()
#> see ?muscData and browseVignettes('muscData') for documentation
#> loading from cache
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
#> dim: 35635 29065 
#> metadata(0):
#> assays(1): counts
#> rownames(35635): MIR1302-10 FAM138A ... MT-ND6 MT-CYB
#> rowData names(2): ENSEMBL SYMBOL
#> colnames(29065): AAACATACAATGCC-1 AAACATACATTTCC-1 ... TTTGCATGGTTTGG-1
#>   TTTGCATGTCTTAC-1
#> colData names(5): ind stim cluster cell multiplets
#> reducedDimNames(1): TSNE
#> mainExpName: NULL
#> altExpNames(0):

counts(sce)[1:4,1:4]
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>              AAACATACAATGCC-1 AAACATACATTTCC-1 AAACATACCAGAAA-1
#> MIR1302-10                  .                .                .
#> FAM138A                     .                .                .
#> OR4F5                       .                .                .
#> RP11-34P13.7                .                .                .
#>              AAACATACCAGCTA-1
#> MIR1302-10                  .
#> FAM138A                     .
#> OR4F5                       .
#> RP11-34P13.7                .

head(colData(sce))
#> DataFrame with 6 rows and 5 columns
#>                        ind     stim   cluster            cell multiplets
#>                  <integer> <factor> <integer>        <factor>   <factor>
#> AAACATACAATGCC-1       107     ctrl         5 CD4 T cells        doublet
#> AAACATACATTTCC-1      1016     ctrl         9 CD14+ Monocytes    singlet
#> AAACATACCAGAAA-1      1256     ctrl         9 CD14+ Monocytes    singlet
#> AAACATACCAGCTA-1      1256     ctrl         9 CD14+ Monocytes    doublet
#> AAACATACCATGCA-1      1488     ctrl         3 CD4 T cells        singlet
#> AAACATACCTCGCT-1      1256     ctrl         9 CD14+ Monocytes    singlet

head(rowData(sce))
#> DataFrame with 6 rows and 2 columns
#>                      ENSEMBL       SYMBOL
#>                  <character>  <character>
#> MIR1302-10   ENSG00000243485   MIR1302-10
#> FAM138A      ENSG00000237613      FAM138A
#> OR4F5        ENSG00000186092        OR4F5
#> RP11-34P13.7 ENSG00000238009 RP11-34P13.7
#> RP11-34P13.8 ENSG00000239945 RP11-34P13.8
#> AL627309.1   ENSG00000237683   AL627309.1
```

For convenience, column data (data describing the cells) can be accessed
using the ‘\$’ operator.

``` r


table(sce$cluster)
#> 
#>    1    2    3    4    5    6    7    8    9   10 
#>  982  115 1141 3232 5755 5179 1905 3569 3921 3266

table(sce$ind,sce$stim)
#>       
#>        ctrl stim
#>   101  1027 1363
#>   107   645  641
#>   1015 3138 2703
#>   1016 2180 1998
#>   1039  552  797
#>   1244 2250 1755
#>   1256 2451 2275
#>   1488 2376 2914
```

We can also store transformations of the raw data, e.g. normalised
counts.

``` r


assays(sce)
#> List of length 1
#> names(1): counts
sce<-logNormCounts(sce)
#> Warning in .library_size_factors(assay(x, assay.type), ...): 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning in .local(x, ...): 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
assays(sce)
#> List of length 2
#> names(2): counts logcounts
```

Aaaand we can also store embeddings of the data, such as from PCA, TSNE
or UMAP.

``` r


reducedDimNames(sce)
#> [1] "TSNE"
```

Finally, we can store alternative experiments, or modalities, which were
measured on the same cells. For example,
[CITE-seq](https://emea.illumina.com/techniques/sequencing/rna-sequencing/cite-seq.html)
data measures gene and surface protein expression simultaneously.

Why might we want to treat these separately rather than appending them
to the main counts?

``` r


altExpNames(sce)
#> character(0)
```

## Exercise

- Can you extract basic summary statistics from the metadata?
- Can you extract the counts data?
- Can you describe the data?

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
#>  [1] muscData_1.27.1             ExperimentHub_3.3.2        
#>  [3] AnnotationHub_4.3.2         BiocFileCache_3.3.0        
#>  [5] dbplyr_2.6.0                scuttle_1.23.2             
#>  [7] SingleCellExperiment_1.35.2 SummarizedExperiment_1.43.0
#>  [9] Biobase_2.73.2              GenomicRanges_1.65.4       
#> [11] Seqinfo_1.3.2               IRanges_2.47.5             
#> [13] S4Vectors_0.51.10           BiocGenerics_0.59.12       
#> [15] generics_0.1.4              MatrixGenerics_1.25.0      
#> [17] matrixStats_1.5.0           BiocStyle_2.41.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1     dplyr_1.2.1          blob_1.3.0          
#>  [4] filelock_1.0.3       Biostrings_2.81.9    fastmap_1.2.0       
#>  [7] digest_0.6.39        lifecycle_1.0.5      KEGGREST_1.53.6     
#> [10] RSQLite_3.53.3       magrittr_2.0.5       compiler_4.6.1      
#> [13] rlang_1.3.0          sass_0.4.10          tools_4.6.1         
#> [16] yaml_2.3.12          knitr_1.52           S4Arrays_1.13.2     
#> [19] htmlwidgets_1.6.4    bit_4.6.0            curl_8.0.0          
#> [22] DelayedArray_0.39.8  abind_1.4-8          BiocParallel_1.47.0 
#> [25] withr_3.0.3          purrr_1.2.2          desc_1.4.3          
#> [28] grid_4.6.1           beachmat_2.29.3      cli_3.6.6           
#> [31] rmarkdown_2.32       crayon_1.5.3         ragg_1.5.2          
#> [34] otel_0.2.0           httr_1.4.9           BiocBaseUtils_1.15.1
#> [37] DBI_1.3.0            cachem_1.1.0         parallel_4.6.1      
#> [40] AnnotationDbi_1.75.2 BiocManager_1.30.27  XVector_0.53.0      
#> [43] vctrs_0.7.3          Matrix_1.7-6         jsonlite_2.0.0      
#> [46] bookdown_0.48        bit64_4.8.6          systemfonts_1.3.2   
#> [49] jquerylib_0.1.4      glue_1.8.1           pkgdown_2.2.1       
#> [52] codetools_0.2-20     BiocVersion_3.24.0   tibble_3.3.1        
#> [55] pillar_1.11.1        rappdirs_0.3.4       htmltools_0.5.9     
#> [58] R6_2.6.1             httr2_1.3.0          textshaping_1.0.5   
#> [61] evaluate_1.0.5       lattice_0.23-1       png_0.1-9           
#> [64] memoise_2.0.1        bslib_0.12.0         Rcpp_1.1.2          
#> [67] SparseArray_1.13.4   xfun_0.61            fs_2.1.0            
#> [70] pkgconfig_2.0.3
```
