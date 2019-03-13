
<!-- README.md is generated from README.Rmd. Please edit that file -->
ZIMHMM
======

The goal of ZIMHMM (Zero Inflated Mixed effects Hidden Mark Model) is to provide a peak caller tailored to detect broad enrichment regions from multiple ChIP-seq experimental replicates. This (under development) package implements `ZIMHMM()` as well as the fixed effects version `ZIHMM()`.

Installation
------------

You can install the released version of ZIMHMM from this repository with:

``` r
devtools::install_github("plbaldoni/ZIMHMM")
```

Example
-------

The package contains an example dataset with H3K36me3 ChIP-seq read counts (with input controls) from chromosome 19 of Huvec cell lines tabulated using 500bp windows (ENCODE Projec). It contains only three functions available for the user: 'ZIMHMM', 'ZIHMM', and 'findpeaks.control()'. Please refer to their documentation for details.

An example of ZIMHMM is shown below.

``` r
# Loading example dataset
data(H3K36me3.Huvec)

# ChIP and Control read counts, as well as the model offset.
ChIP = as.matrix(H3K36me3.Huvec[,c("H3K36me3.Huvec.Rep1","H3K36me3.Huvec.Rep2","H3K36me3.Huvec.Rep3")])
Control = log(as.matrix(H3K36me3.Huvec[,c("Control.Huvec.Rep1","Control.Huvec.Rep2","Control.Huvec.Rep3")])+1)
offset = matrix(0,nrow = nrow(ChIP),ncol = ncol(ChIP))

# Calling peaks
peakcall = ZIMHMM(ChIP = ChIP,Control = Control,offset = offset,random = 'intercept',control = findpeaks.control())
```
