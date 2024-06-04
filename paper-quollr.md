---
title: 'quollr: An R Package for Visualizing 2D Models from Nonlinear Dimension Reductions
  in High Dimensional Space'
description: |
  Non-Linear Dimension Reduction (NLDR) techniques have emerged as powerful tools to visualize high-dimensional data in low-diemnsioanl space. However, their complexity and (hyper)parameter choices may lead to distrustful or misleading results. The R package quollr is developed as a new tool to help to determine which method, which (hyper)parameter choice provide the most accurate representation of $p-D$ data. Clustering data from cardinalR package, is used to illustrate the algorithm and its use within the package.
draft: yes
author:
- name: Jayani P.G. Lakshika
  affiliation: Monash University
  address: Department of Econometrics and Business Statistics, VIC 3800 Australia
  url: https://jayanilakshika.netlify.app/
  orcid_id: 0000-0002-6265-6481
  email: \email{jayani.piyadigamage@monash.edu}
- name: Dianne Cook
  affiliation: Monash University
  address: Department of Econometrics and Business Statistics, VIC 3800 Australia
  url: http://www.dicook.org/
  email: dicook@monash.edu
  orcid_id: 0000-0002-3813-7155
- name: Paul Harrison
  affiliation: Monash University
  address: MGBP, BDInstitute, VIC 3800 Australia
  email: paul.harrison@monash.edu
  orcid_id: 0000-0002-3980-268X
- name: Michael Lydeamore
  affiliation: Monash University
  address: Department of Econometrics and Business Statistics, VIC 3800 Australia
  email: michael.lydeamore@monash.edu
  orcid_id: 0000-0001-6515-827X
- name: Thiyanga S. Talagala
  affiliation: University of Sri Jayewardenepura
  address: Department of Statistics, Gangodawila, Nugegoda 10100 Sri Lanka
  url: https://thiyanga.netlify.app/
  email: ttalagala@sjp.ac.lk
  orcid_id: 0000-0002-0656-9789
type: package
creative_commons: CC BY
date: '2024-06-04'
preamble: |
  \usepackage{amsmath} \usepackage{array}
output:
  distill::distill_article:
    css: style.css
    keep_md: yes
bibliography: RJreferences.bib
editor_options:
  chunk_output_type: inline
journal:
  title: The R Journal
  issn: 2073-4859
  firstpage: 1.0
  lastpage: ~
slug: paper-quollr
pdf_url: paper-quollr.pdf
packages:
  cran:
  - quollr
  - cardinalR
  - testthat
  - knitr
  - rmarkdown
  bioc: []
CTV: ReproducibleResearch
csl: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/rjtools/rjournal.csl

---




<div class="layout-chunk" data-layout="l-body">


</div>


<div class="layout-chunk" data-layout="l-body">


</div>


<!-- 20 pages-->

# Introduction

<!-- research gap: add about hexbin pkg, and emphasize that in our package provide regular hexagons-->
<!-- objective: introduce a new tool to help to determine which method, which parameter choice provide the most useful representation of high-D data.--> 
<!--intro with S-curve with 5 methods-->

This paper presents the R package, `quollr` which introduce a new visual tool in determining which NLDR technique and which (hyper)parameter choice gives most accurate representation of high-dimensional data. The methodology of the algorithm is explained in *cite the methodology paper*. Furthermore, the `quollr` package enables users to perform hexagonal binning [@dan2023], resulting in the generation of regular hexagons. The software is available from the Comprehensive R Archive Network (CARN) at [https://CRAN.R-project.org/package=quollr](https://CRAN.R-project.org/package=quollr).

The paper is organized as follows. In next section, introduces the implementation of `quollr` package on CRAN, including demonstration of the package's key functions and visualization capabilities. We illustrate the algorithm's functionality to study clustering data set in **Application** section, and describe a visual heuristic to describe parameter selection. Finally, we give a brief conclusion of the paper and discuss potential opportunities for use of our algorithm.

# Implementation

The package can be installed from CRAN:

```r
install.packages("quollr")
```

The development version can be installed from GitHub:

```r
devtools::install_github("JayaniLakshika/quollr")
```

## Package dependencies

Understanding the dependencies of the `quollr` package is essential for smooth operation and error prevention. The following dependencies refer to the other R packages that `quollr` relies on to execute its functions effectively. 

<div class="layout-chunk" data-layout="l-body">

```
$quollr
 [1] "dplyr"       "ggplot2"     "grid"        "interp"     
 [5] "langevitour" "proxy"       "rlang"       "rsample"    
 [9] "stats"       "tibble"      "tidyselect" 
```

</div>


## Usage

<!-- add about main function to gen 2D and highD models-->
<!--Discuss the model can be generated with bin centroids or bin means-->

## Constructing the 2-D model

### Scaling the data

### Computing hexagon grid configurations

### Binning the data

### Indicating neighbors by line segments connecting centroids

## Lifting the model into high dimensions

## Model parameters

<!--discuss about default settings-->

## Prediction

## Compute residuals and Mean Square Error (MSE)

## Visualizations

### $2\text{-}D$ model visualization

### $p\text{-}D$ model visualization

## Tests

All functions have tests written and implemented using the \CRANpkg{testthat} [@testthat] in R.

These tests illuminated the issues that allowed us to make meaningful changes and understand some pitfalls of the package.

<!--discuss a test with how a point in a intersection of two hexagonal going to evaluate-->

# Application

<!--with prism data-->

# Discussion

This paper presents the R package `quollr` to develop a way to take the fitted model, as represented by the positions of points in 2D, and turn it into a high-dimensional wireframe to overlay on the data, viewing it with a tour.

The paper includes a clustering example to illustrate how `quollr` is useful to assess which NLDR technique and which (hyper)parameter choice gives the most accurate representation. In addition, how to select parameters for hexagonal binning and fitting model are explained.

Possible future improvements would be...<!--assess the preservation of local and glocal structure w.r.t 2D and high-D distance comparison--> 

This new tool provides an effective start point for automatically creating regular hexagons and help to evaluate which NLDR technique and which hyperparameter choice gives the most accurate representation of $p-D$ data.

# Acknowledgements

This article is created using \CRANpkg{knitr} [@knitr] and \CRANpkg{rmarkdown} [@rmarkdown] in R with the `rjtools::rjournal_article` template. The source code for reproducing this paper can be found at: <https://github.com/JayaniLakshika/paper-quollr>.
```{.r .distill-force-highlighting-css}
```


## CRAN packages used {.appendix}

[quollr](https://cran.r-project.org/package=quollr), [cardinalR](https://cran.r-project.org/package=cardinalR), [testthat](https://cran.r-project.org/package=testthat), [knitr](https://cran.r-project.org/package=knitr), [rmarkdown](https://cran.r-project.org/package=rmarkdown)

## CRAN Task Views implied by cited packages {.appendix}

[ReproducibleResearch](https://cran.r-project.org/view=ReproducibleResearch)




