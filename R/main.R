# This file contains a single function for peak calling with arguments:
# - random: must be either 'none', 'slope', or 'intercept' representing which random effects to be included in the model.
# - offset: either 'sum' or a vector of length equal to the number of experiments with the sequencing depth offset. If a vector of zeroes, no offset is included in the model.
# - group: logical indicating whether to included group effect in the model.
# - quiet: logical indicating whether to output general messages.

#library(Rcpp)
#library(SummarizedExperiment)
#sourceCpp(paste0('/Users/pedroluizbaldoni/Dropbox/PhD/Research/Project1/','Package/Functions/code.cpp'))
#source(paste0('/Users/pedroluizbaldoni/Dropbox/PhD/Research/Project1/','Package/Functions/Base.R'))

findpeaks.control = function(epsilon.em=c(1e-3,1e-3,1e-3,1e-3),epsilon.inner.em=1e-03,maxit.em=500,
                             minit.em=3,gap.em=3,maxcount.em=3,max.phi=1e3,min.sigma2=1e-08,
                             max.sigma2=10,maxcount.inner.em=50,criterion='MULTI',
                             min.zero=.Machine$double.xmin,pcut=0.05,
                             # laptype='none',lapcut=1,lapchr='chr19',chrtb=NULL,
                             quiet=F){
    # if (lapcut<1 & is.null(chrtb)){stop('Provide table with chromosome name, start, stop, and lengths')}
    #if (!(length(epsilon.em)==3) & criterion=='MULTI'){stop("specify criterion for MRCPE, ARCEL, and ACC")}
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
         #laptype=laptype,lapcut=lapcut,lapchr=lapchr,chrtb=chrtb,
         quiet=quiet)
}

findpeaks = function(object,subset,random='none',offset='sum',control=findpeaks.control()){
  ### Checking for input inconsistencies ###
  if(missing(object) | class(object)!='RangedSummarizedExperiment' | !('ChIP'%in%names(object@assays))){stop('Provide valid RangedSummarizedExperiment object.')}
  if(!(subset%in%paste0("chr", c(1:22, "X", "Y")))){stop(paste0('Subset must be a valid chromosome: ',paste0('(',paste0("chr", c(1:22, "X", "Y"),collapse=', '),')')))}
  if(!(random %in% c('none','slope','intercept'))){stop('Provide valid random argument.')}
  if(!(all(offset=='sum') | all(length(offset)==ncol(object) & class(offset)=='numeric'))){stop('Provide valid offset object.')}
  if(!is.list(control) | !all((names(control)%in%names(findpeaks.control())))){stop('Provide valid list of control elements.')}
    if(random %in% c('slope','intercept') & ncol(assay(object,'ChIP'))==1){stop('To fit a random effects model, the number of experiments must be > 1.')}

  ### Assigning objects ###
  control.user = replace(findpeaks.control(),match(names(control),names(findpeaks.control())),control)
  if(!is.null(subset)){object = object[seqnames(object)%in%subset]}

  ChIP = assay(object,'ChIP')
  if('Control' %in% names(object@assays)){
      Control = log(assay(object,'Control')+1)
  } else {
      if(random=='slope'){stop('Provide control experiments to fit a random slope model')}
      Control = NULL
  }
  if(all(offset=='sum')){
      offset = matrix(rep(as.numeric(scale(log(colSums(ChIP)),scale=F)),each=nrow(assay(object))),ncol=ncol(assay(object)))
  } else{
  offset = matrix(rep(offset,each=nrow(assay(object))),ncol=ncol(assay(object)))
  }

  ### Fitting the model ###
  if(random == 'none'){
      model = ZIHMM(ChIP=ChIP,Control=Control,offset=offset,control=control.user)
  } else{
      model = ZIMHMM(ChIP=ChIP,Control=Control,offset=offset,random=random,control=control.user)
  }
  model$Ranges = rowRanges(object)
  model$Ranges$Viterbi = model$Viterbi; model$Ranges$PostProb = model$Prob[,1]

  return(model)
}

# Example
# object
#
# zihmm = findpeaks(object=object,subset='chr22',random='none',offset='sum')
# zimhmm.rs = findpeaks(object=object,subset='chr22',random='slope',offset='sum')
# zimhmm.ri = findpeaks(object=object,subset='chr22',random='intercept',offset=c(0,0,0))
#
# save(zihmm,zimhmm.rs,zimhmm.ri,file='./Vignete_files/Main.RData')
