
<!-- README.md is generated from README.Rmd. Please edit that file -->
ZIMHMM
======

ZIMHMM (Zero Inflated Mixed effects Hidden Mark Model) is a package with a peak caller to detect broad enrichment regions from multiple ChIP-seq experimental replicates. `ZIMHMM()` models the zero-inflation of background counts (ZI), accounts for replicate-specific differences via a mixed effects model (M), and ensures that broad regions of enrichment are detected by fitting a hidden Markov model (HMM). This package also contains `ZIHMM()`, a fixed effects version of the peak caller. The function `findpeaks.control()` provides parameters that control the Expectation-Maximization (EM) algorithm from the presented peak callers. Please refer to their documentation (e.g. `?ZIMHMM::ZIMHMM`) for details.

Installation
------------

You can install the released version of ZIMHMM from this repository with:

``` r
devtools::install_github("plbaldoni/ZIMHMM")
```

Example
-------

The package contains an example dataset with H3K36me3 ChIP-seq (and input control) read counts from chromosome 19 of three technical replicates of Huvec cell line (ENCODE Consortium). Reads were tabulated using 500bp non-overlapping windows. An example of ZIMHMM is shown below.

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
