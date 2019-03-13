#' Control parameters for ZIMHMM and ZIHMM
#'
#' This function passes controlling parameters for ZIMHMM and ZIHMM. Most of these parameters control the EM algorithm
#'
#' @param epsilon.em Either a positive value or a vector of size 4 with the convergence tolerance values for the EM algorithm (see 'criterion' below). Default is c(1e-3,1e-3,1e-3,1e-3)
#' @param epsilon.inner.em positive convergence tolerance for the conditional maximizations used in ZIMHMM
#' @param maxit.em integer giving the maximum number of EM iterations (default 500)
#' @param minit.em integer giving the minimum number of EM iterations to start evaluating the convergence (default 3)
#' @param gap.em integer giving the number of EM iterations apart to compute the convergence criterion (default 3)
#' @param maxcount.em integer giving the number of consecutive EM iterations satisfying the convergence criterion  in order to stop the algorithm (default 3)
#' @param max.phi maximum positive value allowed for the dispersion parameters (default 1000)
#' @param min.sigma2 minimum positive value allowed for the variance component (default 1e-08)
#' @param max.sigma2 maximum positive value allowed for the variance component (default 10)
#' @param maxcount.inner.em integer giving the maxium number of conditional maximizations used in ZIMHMM (default 50)
#' @param criterion convergence criterion: either "MRCPE" (maximum absolute relative change in parameter estimates), "MACPE" (maximum absolute change of parameter estimates),
#' "ARCEL" (absolute relative change of the Q-function), "ACC" (agreement of Viterbi peak calls),
#' or "MULTI" (simultaneously check for MRCPE, MACPE, ARCEL, and ACC).
#' For ACC, it computes the percentage of windows falling in the main diagonal of a 2 by 2 table of Viterbi predictions 'gap.em' iterations apart. Default is "MULTI"
#' @param min.zero minimum positive value allowed in computations to avoid having zeros (default is .Machine$double.xmin)
#' @param pcut cutoff for rejection controlled EM algorithm (default 0.05)
#' @param quiet whether to print messages (default F)
#'
#' @return A list with components equal to the arguments
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/ZIMHMM}
#'
#' @examples
#' data(H3K36me3.Huvec)
#' ChIP = as.matrix(H3K36me3.Huvec[,c("H3K36me3.Huvec.Rep1","H3K36me3.Huvec.Rep2","H3K36me3.Huvec.Rep3")])
#' Control = log(as.matrix(H3K36me3.Huvec[,c("Control.Huvec.Rep1","Control.Huvec.Rep2","Control.Huvec.Rep3")])+1)
#' offset = matrix(0,nrow = nrow(ChIP),ncol = ncol(ChIP))
#' # Setting maxit.em = 100 (no more than 100 EM iterations)
#' control = findpeaks.control(maxit.em = 100)
#' ZIMHMM(ChIP = ChIP,Control = Control,offset = offset,random = 'intercept',control = control)
#'
#' @useDynLib ZIMHMM
#' @importFrom Rcpp evalCpp
#' @export
#'
findpeaks.control = function(epsilon.em=c(1e-3,1e-3,1e-3,1e-3),epsilon.inner.em=1e-03,maxit.em=500,
                             minit.em=3,gap.em=3,maxcount.em=3,max.phi=1e3,min.sigma2=1e-08,
                             max.sigma2=10,maxcount.inner.em=50,criterion='MULTI',
                             min.zero=.Machine$double.xmin,pcut=0.05,
                             quiet=F){
    if (!is.numeric(epsilon.em) || epsilon.em <= 0){stop("value of 'epsilon.em' must be > 0")}
    if (!maxit.em%%1==0 || maxit.em <= 0){stop("value of 'maxit.em' must be a positive integer")}
    if (!minit.em%%1==0 || minit.em <= 0){stop("value of 'minit.em' must be a positive integer")}
    if (!gap.em%%1==0 || gap.em <= 0 || gap.em>minit.em){stop("value of 'gap.em' must be a positive integer <= minit.em")}
    if (!maxcount.em%%1==0 || maxcount.em <= 0){stop("value of 'maxcount.em' must be a positive integer")}
    if (!is.numeric(max.phi) || max.phi <= 0){stop("value of 'max.phi' must be > 0")}
    if (!is.numeric(max.sigma2) || max.sigma2 <= 0){stop("value of 'max.sigma2' must be > 0")}
    if (!maxcount.inner.em%%1==0 || maxcount.inner.em <= 0){stop("value of 'maxcount.inner.em' must be a positive integer")}
    if (!is.logical(quiet)){stop("value of 'quiet' must be logical TRUE or FALSE")}
    list(epsilon.em=epsilon.em,epsilon.inner.em=epsilon.inner.em,maxit.em=maxit.em,minit.em=minit.em,
         gap.em=gap.em,maxcount.em=maxcount.em,max.phi=max.phi,min.sigma2=min.sigma2,
         max.sigma2=max.sigma2,maxcount.inner.em=maxcount.inner.em,criterion=criterion,
         min.zero=min.zero,pcut=pcut,
         quiet=quiet)
}
