# library(SummarizedExperiment)
# library(chipseq)
# library(Rsamtools)
# library(rtracklayer)
# library(bamsignals)
# library(csaw)
# library(ChIPQC)
#
# # Loading genome
# setwd('/Users/pedroluizbaldoni/Dropbox/PhD/Research/Project1/Package/')
# load('./Misc/hg19.RData')

# List of functions
# - readcounts: this function takes as argument a set of ChIP bam files, a set of control bam files asssociated
#     with the ChIP experiments, a vector of condition/replicate labels and an integer with the bin size. It outputs a
#     SummarizedExperiment element with assays containing the tabulated read counts.

readcounts = function(sample,control=NULL,condition,replicate,size){
  if(is.null(control)){if(!all(sapply(list(sample,condition,replicate),length)==length(list(sample,condition,replicate)[[1]]))){stop('Error: different input lengths.')}}
  if(!is.null(control)){if(all(sapply(list(sample,control,condition,replicate),length)==length(list(sample,control,condition,replicate)[[1]]))==F){stop('Error: different input lengths.')}}
  if(missing(size)){stop('Error: provide valid window size.')}

  N = length(sample)

  cat('Estimating fragment length...\n')
  flength = unlist(lapply(lapply(sample,csaw::correlateReads,param=csaw::readParam(discard=hg19.discard)),csaw::maximizeCcf))
  cat('Tiling the genome...\n')
  gr.tile = with(data.frame(GenomicRanges::tile(hg19,width=size))[,c('seqnames','start','end')],GenomicRanges::GRanges(seqnames,IRanges::IRanges(start,end)))

  ChIP.counts = matrix(NA,nrow=length(gr.tile),ncol=N)
  if(!is.null(control)){Control.counts = matrix(NA,nrow=length(gr.tile),ncol=N)}
  for(i in 1:N){cat(paste0('\rTabulating reads of ',condition[i],', replicate number ',replicate[i],'...\n'))
    #Check for valid chromosomes in Bamfiles
      # The following two lines cause error because of seqlevels/seqnames:
      #Error in match(x, table, nomatch = 0L) :
      #'match' requires vector arguments
    valid.chip = as.logical(GenomeInfoDb::seqnames(gr.tile) %in% GenomeInfoDb::seqlevels(Rsamtools::BamFile(sample[i])))
    #if(!is.null(control)){valid.control = as.logical(GenomeInfoDb::seqnames(gr.tile) %in% GenomeInfoDb::seqlevels(Rsamtools::BamFile(control[i])))}

    #Tabulating
    ChIP.counts[valid.chip,i] <- (bamsignals::bamCount(bampath=sample[i],gr.tile[valid.chip],verbose=F,shift=flength[i]/2))
    if(!is.null(control)){Control.counts[valid.control,i] <- (bamsignals::bamCount(bampath=control[i],gr.tile[valid.control],verbose=F))}
  }
  if(is.null(control)){
    outp = SummarizedExperiment::SummarizedExperiment(assays=list(ChIP=ChIP.counts),rowRanges=gr.tile,colData=data.frame(Sample=1:N,Condition=condition,Replicate=replicate,Fragment=flength))
  } else {
    outp = SummarizedExperiment::SummarizedExperiment(assays=list(ChIP=ChIP.counts,Control=Control.counts),rowRanges=gr.tile,colData=data.frame(Sample=1:N,Condition=condition,Replicate=replicate,Fragment=flength))
  }
  cat(paste0('Done!'))
  return(outp)
}

# Example
# datawd = '/Users/pedroluizbaldoni/Dropbox/PhD/Research/Project1/Package/Data/'
#
# sample = list.files(paste0(datawd,'ChIP/'),pattern='*.bam$',full.names=T)
# control = list.files(paste0(datawd,'Control/'),pattern='*.bam$',full.names=T)
# condition = c('Huvec','Huvec','Huvec')
# replicate = c(1,2,3)
# size=500
#
# object = readcounts(sample=sample,control=control,condition=condition,replicate=replicate,size=size)
