<!-- README.md is generated from README.Rmd. Please edit that file -->
rrelrates
=========

[![Travis-CI Build Status](https://travis-ci.org/rachelss/rrelrates.png?branch=master)](https://travis-ci.org/rachelss/rrelrates)

Produce relative rates on from trees using the reltime method of Tamura et al. 2012.

This package is a work in progress. Currently it is limited to primates and data available in ensembl. You can install the package via `devtools`

``` r
devtools::install_github("rachelss/rrelrates")
```

Examples
--------

### Mean and Variance of the relative rates

``` r
library(rrelrates)
#> Loading required package: rensembl
#> Loading required package: phytools
#> Loading required package: ape
#> Loading required package: maps
#> Loading required package: phangorn
#> Loading required package: plotrix
geneinfo("CCND1")
#> [1] "CCND1"
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

    #> $CCND1_ENSPPYP00000003434_Pabe
    #> [1] 2.232375
    #> 
    #> $mean
    #> [1] 0.9768734
    #> 
    #> $var
    #> [1] 0.3532899
