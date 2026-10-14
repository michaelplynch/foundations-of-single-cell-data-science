# Foundations of single-cell data science

This workshop introduces the basics of single-cell RNA sequencing (gene
expression) and a minimal analysis workflow to an audience with a
background in data science, machine learning, and statistics. The
workshop is delivered as an installable R package, with a rendered [html
website](https://michaelplynch.github.io/foundations-of-single-cell-data-science/),
and Docker container. An emphasis will be placed on the how and why, not
just the what, of each step.

## Broad Goals

- Understand the structure of single-cell RNA-seq data.
- Access the Bioconductor ecosystem.
- Perform basic single-cell analysis.
- Understand factors impacting analysis reliability.

## This is

- a crash course in the key steps of a single-cell RNA-seq analysis.
- a discussion on how to approach new problems generally in DS.
- considerations for effective analyses.

## This isn’t

- a deep dive of single-cell biology/immunology.
- a demonstration of the only/best way to complete a single-cell
  analysis.

## Workshop breakdown

| Activity                    | Notebook                    | Time          |
|-----------------------------|-----------------------------|---------------|
| Setting up the environment  |                             | 10m           |
| Background and planning     | a_background.Rmd            | 15m           |
| Bioc packages and object    | b_ecosystem_objects.Rmd     | 15m           |
| Single-sample workflow      | c_singlesample_workflow.Rmd | 60m           |
| Reflection & Best Practices | c_singlesample_workflow.Rmd | 15m           |
| Multi-sample workflow       | d_multisample_workflow.Rmd  | self-directed |

## Installation

Rendered pkgdown website available here:
<https://michaelplynch.github.io/foundations-of-single-cell-data-science/>

The workshop can be installed locally using the below commands, from a
fresh R session.

``` r


## Check R version
if (getRversion() < "4.6.0") stop("This workshop needs R 4.6.0 or later (you have ", getRversion(), "). Please update R and re-run.")

## install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

options(repos = BiocManager::repositories(version="3.23"))

## install remotes
install.packages(c('remotes','usethis'))

## install workshop package (takes care of dependencies)
remotes::install_github("michaelplynch/foundations-of-single-cell-data-science",dependencies=TRUE, upgrade='always')

## Course notebooks (same release); answer the prompts, notebooks are in vignettes/
usethis::use_course(
  "https://github.com/michaelplynch/foundations-of-single-cell-data-science/archive/refs/tags/v1.0.0.zip",
  destdir = "~"
)

## Pre-download the dataset and check everything loads
sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
sce <- muscData::Kang18_8vs8()
pkgs <- c("SingleCellExperiment", "scater", "scrapper", "dittoSeq", "BiocStyle","UCell","mclust")
invisible(lapply(pkgs, library, character.only = TRUE))
BiocManager::valid()
message("Setup complete - see you at the workshop!")
```

Tested on Bioconductor 3.23 and R Version 4.6.x.

### Running through Docker

## Running with Docker (alternative to a local install)

Instead of installing packages locally, you can use the pre-built
container. You’ll need [Docker
Desktop](https://www.docker.com/products/docker-desktop/) installed and
running.

``` bash
docker run --rm -p 8787:8787 -e PASSWORD=bioc ghcr.io/michaelplynch/foundations-of-single-cell-data-science:latest
```

Then open <http://localhost:8787> in your browser and log in with
username `rstudio` and password `bioc`. The notebooks are in the
`vignettes/` folder.

- **Apple Silicon Macs (M1–M4):** add `--platform linux/amd64` after
  `docker run`.
- **The first run downloads several GB,** so do it before the workshop.
- **Work is lost when the container stops.**

## Pre-requisites

- Basic knowledge of R syntax
- Basic machine learning and statistics
- Previous data science/analysis experience

## Useful resources

R

- [Orchestrating Single-cell Analysis
  (OSCA)](https://bioconductor.org/books/release/OSCA/)
- [Orchestrating Spatial Transcriptomic Analysis
  (OSTA)](https://bioconductor.org/books/release/OSTA/), a related and
  upcoming research area.

Python

- [Single-cell best practices](https://www.sc-best-practices.org/)

## Participation

Interactive discussion on analysis strategy and results, follow-along
coding, reflection on analysis strengths and limitations.

## Data

Data used are from the *4k PBMCs from a Healthy Donor* dataset released
by 10x Genomics (2017), distributed via the TENxPBMCData Bioconductor
package, and from Kang HM, et al. (2018). Multiplexed droplet
single-cell RNA-sequencing using natural genetic variation. *Nature
Biotechnology* 36:89–94.
[doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042). PBMCs from 8
lupus patients, control and IFN-β stimulated, distributed via the
muscData Bioconductor package.
