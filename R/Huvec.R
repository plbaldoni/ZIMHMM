#' ChIP-seq read counts of H3K36me3 (and input control) from the Huvec cell line of
#' chromosome 19.
#'
#' Data from ENCODE project of H3K36me3 (and input control) ChIP-seq experiments (3 replicates each) from Huvec cells of chromosome 19.
#' PCR duplicates were removed with 'samtools markdup'. Low quality reads were removed with
#' 'samtools view -q 10'. The fragment length was estimated using correlateReads
#' and maximizeCcf (from package csaw) after excluding blacklisted positions.
#' Read counts were tabulated using bamCount (from package bamsignals).
#'
#' @docType data
#'
#' @usage data(Huvec)
#'
#' @format RangedSummarizedExperiment
#'
#' @keywords datasets
#'
#' @references ENCODE Project Consortium, 2012. An integrated encyclopedia of DNA elements in the human genome. Nature, 489(7414), p.57.
#'
#' @source
#' https://www.ncbi.nlm.nih.gov/gds/?term=GSM733757
#' https://www.ncbi.nlm.nih.gov/gds/?term=GSM733715
#'
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' data(Huvec)
#' ChIP = SummarizedExperiment::assay(Huvec,'ChIP')
"Huvec"
