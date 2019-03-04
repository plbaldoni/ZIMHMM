# library(ggplot2)
# library(reshape)
# library(GenomicRanges)
# library(SummarizedExperiment)

calls.df = function(gr.peaks,gr.background,ylim,xlim,method,chromosome,tracknum){
    peaks=intersect(GRanges(chromosome,IRanges(xlim[1],xlim[2])),gr.peaks)
    background=intersect(GRanges(chromosome,IRanges(xlim[1],xlim[2])),gr.background)

    return(rbind(data.frame(x1 = start(peaks),y1 = ylim[1],x2 = end(peaks),y2 = ylim[1],Method='Enrichment',variable=tracknum),
                 data.frame(x1 = start(background),y1 = ylim[1],x2 = end(background),y2 = ylim[1],Method='Background',variable=tracknum)))
}

plotpeaks = function(name,peaks,object,subset=NULL,zoom=50,pct.ylim=0.10,gap=10,range=NULL,ranges=NULL){
    pct.gr = zoom/100
    if(!is.null(subset)){object = object[seqnames(object)%in%subset]}
    #Selects start/stop position of the peak
    start.gr = start(peaks[peaks$name==name])
    stop.gr = end(peaks[peaks$name==name])
    #Selects the relative start/stop position of the data
    start.dat =  start(rowRanges(object)[which.min(abs(start(rowRanges(object))-start.gr))])
    stop.dat = end(rowRanges(object)[which.min(abs(end(rowRanges(object))-stop.gr))])
    width.dat = stop.dat - start.dat
    #Enlarging the start/stop position of the data
    startpct.dat = start.dat-width.dat*pct.gr
    stoppct.dat = stop.dat+width.dat*pct.gr
    #Subsetting read counts
    subdat = object[which.min(abs(start(rowRanges(object))-startpct.dat)):which.min(abs(end(rowRanges(object))-stoppct.dat))]
    subdatdf = data.frame(chr=seqnames(subdat),start=start(subdat),stop=end(subdat),y=assay(subdat))
    subdat.melt = melt(subdatdf[,c('chr','start','stop',paste0('y.',1:ncol(assay(subdat))))],id=c('chr','start','stop'))
    levels(subdat.melt$variable) = paste(colData(object)$Condition,colData(object)$Replicate)
    #Defining coordinates of peaks
    xlim1 = subdatdf$start[1]
    xlim2 = subdatdf$stop[nrow(subdatdf)]
    #Defining peaks now
    peaks.track = calls.df(gr.peaks=peaks,gr.background=setdiff(rowRanges(object),peaks),ylim=max(subdat.melt$value)*(1+pct.ylim),xlim=c(xlim1,xlim2),chromosome=as.character(runValue(seqnames(object))),
                           tracknum = levels(subdat.melt$variable)[1])
    #Creating plot now
    p = ggplot(subdat.melt, aes(x = start, y = value)) + geom_line(color='black') + facet_grid(variable ~ .) +
        scale_x_continuous(limits = c(xlim1, xlim2),breaks=c(xlim1, xlim2))+
        theme_bw()+
        theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8),strip.text.y = element_text(size=8))+
        labs(x=paste('Window index on',as.character(runValue(seqnames(object)))),y='Window Read Count')+ theme(legend.position="bottom",legend.direction='horizontal',legend.title = NULL)
    p = p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2,group=variable,colour=Method),size=3, data = peaks.track,inherit.aes=FALSE)+
        scale_color_manual(values = RColorBrewer::brewer.pal(9, "Set1")[c(1,9)])
    p = p +theme(legend.position="bottom", legend.box = "horizontal",legend.title=element_blank())#theme(legend.position="none")
    print(p)
    return(p)
}
