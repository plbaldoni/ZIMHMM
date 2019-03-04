# library(rtracklayer)

getalpha = function(alpha,p,fdr){sum(p*(p<=alpha))/sum(p<=alpha)-fdr}

summarizepeaks = function(peaks,file,trackname='User Track',method='Viterbi',fun='mean',fdr=0.05)
{
  if(!is.character(file)){stop('The argument "name" must be a character.')}
  if(!method%in%c('Viterbi','FDR')){stop('Method must be either "Viterbi" or "FDR".')}
  if(method=='FDR'){if(!(fdr>0 & fdr<1)){stop('FDR level must be within the (0,1) interval.')}}
  #Reduce peaks
  if(method=='Viterbi'){out = reduce(peaks$Ranges[peaks$Ranges$Viterbi==1])}
  if(method=='FDR'){out = reduce(peaks$Ranges[with(peaks$Ranges,(PostProb<=uniroot(getalpha,interval=c(min(PostProb),max(PostProb)),p=PostProb,fdr=fdr)$root))])}
  #Summarize PostProb
  overlaps = findOverlaps(peaks$Ranges,out)
  PostProb = peaks$Ranges[queryHits(overlaps)]$PostProb
  averagedSignal <- aggregate(PostProb, list(subjectHits(overlaps)),fun)
  out$signalValue = -log10(averagedSignal$x)
  #Organizing output
  out$score = with(out,(1000/(max(signalValue)-min(signalValue)))*(signalValue-max(signalValue))+1000)
  out$pValue = -1
  out$qValue = -1
  #Sorting output
  out = out[order(-out$signalValue)]
  out$name = paste0('Peak',1:length(out))
  #Exporting
  bed = data.frame(chrom=seqnames(out),chromStart=start(out),chromEnd=end(out),name=out$name,score=out$score,strand='.',
                   signalValue=out$signalValue,pValue=out$pValue,qValue=out$qValue)
  write.table(bed,file="temp.bed",row.names=F,col.names=F,quote=F,sep="\t")

  header1 = paste0('track name="ZIMHMM Peaks - ',trackname,'" description="" useScore=1 type=broadPeak description="" useScore=1')
  header2 = paste0('browser position ',bed[1,'chrom'],':',bed[1,'chromStart'],'-',bed[1,'chromEnd'])

  system(paste0("echo '",header2,"' | cat - temp.bed > temp1.bed"))
  system(paste0("echo '",header1,"' | cat - temp1.bed > ",file,'.broadPeak'))
  system('rm temp.bed temp1.bed')

  cat(paste0('A .broadPeak file with peak calls was saved under name: ',file,'.broadPeak'))
  return(out)
}
