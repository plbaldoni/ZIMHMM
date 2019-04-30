#' Create Offsets for Peak Calling
#'
#' This function create offsets for peak calling. The current implementation can output three different set of offsets:
#' total sum of read counts (method="sum"), smoothed MA trend (method="loess", similar to csaw), or median ratio (method="ratio", similar to DESeq2).
#'
#' @param ChIP M*N matrix of ChIP read counts, where M is the number of windows in the analyzed genome and N is the number of experiments and replicates
#' @param method either total sum of read counts (method="sum"), smoothed MA trend (method="loess"), or median ratio (method="ratio").
#' @param span proportion of data to be used in the smoothing (see loessFit from limma, default 0.3)
#'
#' @return M*N matrix of offsets
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/ZIMHMM}
#'
#' @examples
#' data(Huvec)
#' ChIP = SummarizedExperiment::assay(Huvec,'ChIP')
#' offset = createOffset(ChIP,method="loess")
#'
#' @importFrom limma loessFit
#' @export

createOffset = function(ChIP,method='sum',span=0.3){
    N = ncol(ChIP)
    M = nrow(ChIP)
    offset = matrix(NA,nrow=M,ncol=N)
    if(method=='sum'){
        offset = matrix(log(colSums(ChIP+1)),nrow=M,ncol=N,byrow=T)
        return(offset)
    }
    if(method=='loess'){
        log.ChIP = log(ChIP+1)
        avg.ChIP = exp(rowMeans(log.ChIP))
        log.avg.ChIP = log(avg.ChIP)

        log.ChIP.M = log.ChIP - matrix(log.avg.ChIP,nrow=M,ncol=N,byrow = F) #M, from MA plot
        log.ChIP.A = 0.5*(log.ChIP + matrix(log.avg.ChIP,nrow=M,ncol=N,byrow = F)) #A, from MA plot

        # Calculating Loess
        offset = sapply(1:N,function(i){limma::loessFit(y=log.ChIP.M[,i],x=log.ChIP.A[,i],span = span)$fitted})

        return(offset)
    }
    if(method=='ratio'){
        log.ChIP = log(ChIP+1)
        avg.ChIP = exp(rowMeans(log.ChIP))
        log.avg.ChIP = log(avg.ChIP)

        log.ChIP.M = log.ChIP - matrix(log.avg.ChIP,nrow=M,ncol=N,byrow = F) #M, from MA plot
        log.ChIP.A = 0.5*(log.ChIP + matrix(log.avg.ChIP,nrow=M,ncol=N,byrow = F)) #A, from MA plot

        # Calculating offset
        medianRatio = apply(log.ChIP.M,2,median)
        offset = matrix(medianRatio,nrow = M,ncol = N,byrow = T)

        return(offset)
    }
}
