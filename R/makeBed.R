#' Create a BED file from an output of ZIMHMM
#'
#' `makeBed()` imports the output from `ZIMHMM()` and a RangedSummarizedExperiment object dataset to
#' create a BED file to be used in downstream analyses.
#'
#' @param output An output from `ZIMHMM()`
#' @param data A RangedSummarizedExperiment object with the same assay used in `ZIMHMM()`
#' @param file The full path and name of the output file (e.g. './example.bed')
#' @param tracklabel Optional label for the track (default is 'User Track')
#'
#' @return A GRanges object with consensus peaks.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/mixHMM}
#'
#' @examples
#' data(Huvec)
#' ChIP = SummarizedExperiment::assay(Huvec,'ChIP')
#' Control = log(SummarizedExperiment::assay(Huvec,'Control')+1)
#' offset = matrix(0,nrow = nrow(ChIP),ncol = ncol(ChIP))
#' \dontrun{output = ZIMHMM(ChIP = ChIP,Control = Control,offset = offset,random = 'intercept',control = controlPeaks())}
#' \dontrun{makeBed(output,data = Huvec,file = './H3K36me3.bed')}
#'
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom SummarizedExperiment rowRanges width seqnames
#'
#' @export
#'
makeBed = function(output,data,file,tracklabel='User Track')
{
    if(!is.character(file)){stop('The argument "file" must be character.')}
    # Reduce peaks
    out = GenomicRanges::reduce(SummarizedExperiment::rowRanges(data)[output$Viterbi==1])

    #Summarize PostProb
    overlaps = GenomicRanges::findOverlaps(SummarizedExperiment::rowRanges(data),out)
    PostProb = output$Prob$PostProb2[S4Vectors::queryHits(overlaps)]
    averagedSignal <- aggregate(PostProb, list(S4Vectors::subjectHits(overlaps)),'mean')
    out$signalValue = -log10(averagedSignal$x)

    #Organizing output
    out$score = ((1000-500)/(max(out$signalValue)-min(out$signalValue)))*(out$signalValue-max(out$signalValue))+1000
    out$pValue = -1
    out$qValue = -1

    #Sorting output
    out = out[order(SummarizedExperiment::width(out),-out$signalValue)]
    out$name = paste0('Peak',1:length(out))

    #Exporting
    bed = data.frame(chrom=SummarizedExperiment::seqnames(out),chromStart=SummarizedExperiment::start(out),chromEnd=SummarizedExperiment::end(out),name=out$name,score=out$score,strand='.',
                     signalValue=out$signalValue,pValue=out$pValue,qValue=out$qValue)
    write.table(bed,file="./temp.bed",row.names=F,col.names=F,quote=F,sep="\t")

    header1 = paste0('track name="ZIMHMM Peaks - ',tracklabel,'" description="" type=broadPeak useScore=1')
    header2 = paste0('browser position ',bed[1,'chrom'],':',bed[1,'chromStart'],'-',bed[1,'chromEnd'])

    system(paste0("echo '",header2,"' | cat - ./temp.bed > ./temp1.bed"))
    system(paste0("echo '",header1,"' | cat - ./temp1.bed > ",file))
    system('rm ./temp.bed ./temp1.bed')

    cat(paste0('A .bed file with peak calls was saved under name: ',file,'.\n'))
    return(out)
}
