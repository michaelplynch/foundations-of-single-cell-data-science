# Foundations of single-cell data science

This package contains the workshop material for the Foundations of Data Science CRT Single-cell and Bioconductor workshop.
An emphasis will be placed on the how and why, not just the what, of each step.

## This is

- a crash course in the key steps of a single-cell RNA-seq analyses.
- a discussion on how to approach new problems generally in DS.
- considerations for effective analyses.

## This isn't

- a deep dive of single-cell biology/immunology.
- a demonstration of the only/best way to complete a single-cell analysis.

## Workshop breakdown

| Activity                     | Time |
|------------------------------|------|
| Setting up the environment   | 10m  |
| Planning an analysis         | 15m  |
| Bioconductor ecosystem       | 15m  |
| Single-sample workflow       | 60m  |
| Reflection & Best Practices  | 15m  |

## Installation:

Rendered pkgdown website available here: https://michaelplynch.github.io/foundations-of-single-cell-data-science/

Package note

```r
## install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

options(repos = BiocManager::repositories())

## install remotes
install.packages('remotes')

## install workshop package (takes care of dependencies)
remotes::install_github("michaelplynch/foundations-of-single-cell-data-science",dependencies=T)

## download course
usethis::use_course("michaelplynch/foundations-of-single-cell-data-science")

```


Tested on R Version 4.6.1.


### Description

This workshop introduces the basics of single-cell RNA sequencing (gene expression) and a minimal analysis workflow to students with a background in data science, machine learning, and statistics.

### Pre-requisites

* Basic knowledge of R syntax
* Basic machine learning and statistics
* Some previous data science/analysis experience

### Useful resources

- https://bioconductor.org/books/release/OSCA/
- https://bioconductor.org/books/release/OSTA/

### Participation

Interactive discussion, follow-along coding, reflection.




## Below to be moved/tidied up pre workshop
## Map:

### How might a data scientist approach this problem?

1. Understand how the data was generated.
2. Understand the type of analysis (exploratory, hypothesis generation, vs. hypothesis testing).
3. Check the data structure, quality.
4. Decide on a clear and defensible modelling strategy.
5. Validate outputs. Big difference between 'an output' and 'a correct output'. Use domain expertise, independent datasets, stability and reproducibility checks.
6. Communicate clearly. What claims can you make, and what claims can't you make, based on the data?

### Foundations

- Statistical thinking
  - Data distributions
  - Independence
  - Multiple testing
  - Assumptions
- Computational thinking
  - Scaling
  - Seeds, reproducibility
  - Algorithms
- Communication, visualisation
  - Model assumptions
  - Uncertainty
  - Informative vs. misleading
- Tooling
  - Choice of ecosystem
  - Version control
  - Documentation
  - Dependency management (note to check number of packages used).
- Ethics, integrity
  - Reproducibility (again).
  - Choice of reference.
  - Patient privacy and consent.
