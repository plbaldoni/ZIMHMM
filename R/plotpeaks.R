#' Create a plot of ChIP read counts with called peaks and posterior probabilities
#'
#' `plotPeaks()` imports the output from `ZIMHMM()` or `ZIHMM()` and plot the peak
#' calls of a given genomic region.
#'
#' @param output Either an output from `ZIMHMM()` or `ZIHMM()`
#' @param chr A string with the chromosome name to be displayed (such as 'chr19')
#' @param ranges A vector with two elements that specifies the genomic windows to be plotted.
#' @param data A RangedSummarizedExperiment object with the same assay (with the `ChIP` counts) used in `ZIMHMM()`
#' @param peaks An optional vector of length M with zeros (background) and ones (peaks) representing the enriched windows. If not provided,
#' `plotPeaks()` will use the Viterbi sequence from `output`.
#'
#' @return A ggplot
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/ZIMHMM}
#'
#' @examples
#' data(Huvec)
#' ChIP = SummarizedExperiment::assay(Huvec,'ChIP')
#' Control = log(SummarizedExperiment::assay(Huvec,'Control')+1)
#' offset = matrix(log(colSums(ChIP)),nrow = nrow(ChIP),ncol = ncol(ChIP),byrow = TRUE)
#' \dontrun{output = ZIHMM(ChIP = ChIP,Control = Control,offset = offset,control = controlPeaks(epsilon.em = 1e-3,criterion = 'MRCPE'))}
#' \dontrun{plotPeaks(output = output,chr = 'chr19',ranges = c(16898093,17040483),data = Huvec)}
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
#' @importFrom dplyr as.tbl
#' @importFrom plyr mapvalues
#' @importFrom ggpubr ggarrange
#' @importFrom scales comma
#' @importFrom GenomicRanges seqnames
#'
#' @export
#'
plotPeaks = function(output,chr,ranges,data,peaks = NULL){
    if(length(ranges)!=2 | ranges[2]<ranges[1]){stop('Argument ranges must be a vector of two integers (x,y) with y>x')}
    if(!any(GenomicRanges::seqnames(data)==chr)){stop('Chromosome not found')}

    oldranges = ranges
    data <- data[GenomicRanges::seqnames(data)==chr]
    ranges = c(which.min(abs(start(rowRanges(data))-ranges[1])),which.min(abs(end(rowRanges(data))-ranges[2])))
    ranges = ranges[1]:ranges[2]

    datamelt <- data.table::melt(as.data.table(assay(data,'ChIP')),measure.vars=seq_len(ncol(assay(data,'ChIP'))),value.name = 'Counts',variable.name = 'Replicates')
    setattr(datamelt$Replicates,"levels",paste0('Replicate ',1:ncol(ChIP)))
    datamelt[,Window := rep(seq_len(nrow(assay(data,'ChIP'))),ncol(assay(data,'ChIP')))]

    datapeak = data.table::data.table(Window = seq_len(nrow(assay(data,'ChIP'))),Replicates = unique(datamelt$Replicates)[1])
    if(is.null(peaks)){
        datapeak[,Peaks := output$Viterbi]
    } else{
        datapeak[,Peaks := peaks]
    }
    datapeak[Peaks==0,Peaks := NA]

    ### Figure 1: Read Counts and Peak Calls ###
    maxy = max(datamelt[(Window%in%ranges),Counts])*1.1
    color = RColorBrewer::brewer.pal(3,'Set1')[2]
    fig.ChIP = ggplot2::ggplot(data=datamelt[(Window%in%ranges),],ggplot2::aes(x=Window,y=Counts))+
        ggplot2::geom_line()+
        ggplot2::facet_grid(rows=ggplot2::vars(Replicates))+
        ggplot2::geom_segment(inherit.aes=FALSE,data = datapeak[(Window%in%ranges),],ggplot2::aes(x=Window,xend=Window+1,y=maxy*Peaks,yend=maxy*Peaks),size=2,color=color)+
        ggplot2::theme_bw()+
        ggplot2::xlab('Genomic Window')+
        ggplot2::ylab('Read Counts')+
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank())

    ### Figure 2: Post. Probabilities ###
    PostProb = data.table::data.table(Window=start(rowRanges(data)),output$Prob,Label='PP')
    PostProb <- tidyr::gather(dplyr::as.tbl(PostProb),Component,PP,PostProb2)
    PostProb$Component = plyr::mapvalues(PostProb$Component,from=c('PostProb2'),to=c('Enrichment'))

    fig.Prob = ggplot2::ggplot(data=PostProb[PostProb$Window%in%start(rowRanges(data)[ranges]),],ggplot2::aes(x=Window,y=PP,fill=Component))+
        ggplot2::facet_grid(rows=ggplot2::vars(Label))+
        ggplot2::geom_area(position='identity',alpha=1,color=color,fill=color)+
        ggplot2::scale_y_continuous(limits = c(0,1),breaks=c(0,0.5,1))+
        ggplot2::scale_x_continuous(limits = range(PostProb[PostProb$Window%in%start(rowRanges(data)[ranges]),'Window']),
                                    labels = scales::comma)+
        ggplot2::xlab(paste0('Genomic Window'))+
        ggplot2::ylab('Post. Prob.')+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = 'bottom',legend.direction = 'horizontal',
                       strip.text.y = ggplot2::element_text(colour = ggplot2::alpha('grey',0.0)))+
        ggplot2::theme(legend.position="none")

    suppressWarnings({fig <- ggpubr::ggarrange(fig.ChIP,fig.Prob,ncol=1,nrow=2,heights = c(0.8,0.2))})

    return(fig)
}
