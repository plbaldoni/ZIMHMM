#' Hidden Markov Model (HMM)
#'
#' This function runs a 2-state HMM with negative binomial emissions. It is used for parameter initialization in ZIMHMM and ZIHMM.
#'
#' @param ChIP.init M*N matrix of ChIP read counts, where M is the number of windows in the analyzed genome and N is the number of replicates
#' @param Control.init M*N matrix of log-transformed Control read counts
#' @param offset.init M*N matrix of offsets. If no offset is used, use offset = matrix(0,nrow=M,ncol=N)
#' @param pcut cutoff for rejection controlled EM algorithm (default is 0.05)
#' @param epsilon.em A positive value with the convergence tolerance value for the EM algorithm (default is 1e-3)
#' @param maxit.em integer giving the maximum number of EM iterations (default 5)
#' @param minit.em integer giving the minimum number of EM iterations to start evaluating the convergence (default 3)
#' @param gap.em integer giving the number of EM iterations apart to compute the convergence criterion (default 3)
#' @param maxcount.em integer giving the number of consecutive EM iterations satisfying the convergence criterion  in order to stop the algorithm (default 3)
#' @param max.phi maximum positive value allowed for the dispersion parameters (default 1000)
#' @param min.zero minimum positive value allowed in computations to avoid having zeros (default is .Machine$double.xmin)
#' @param quant quantile of the distribution of the score of scaled ChIP counts to define the very first set of enrichment and background windows to begin the EM algorithm (default is 0.75).
#' Scaled windows with score below the quant-tile are defined as background to start the EM algorithm.
#' @param quiet whether to print messages (default F)
#'
#' @return A list with components
#' \item{Pi}{Vector of initial probabilities of the HMM}
#' \item{Gamma}{Matrix of transition probabilities of the HMM}
#' \item{Psi}{Vector of component-specific parameters of the HMM}
#' \item{Prob}{Mx2 Matrix with posterior probabilities}
#' \item{LogF}{Mx2 Matrix with log-forward probabilities}
#' \item{LogB}{Mx2 Matrix with log-backward probabilities}
#' \item{Loglik}{Mx2 Matrix with window-based probabilities}
#' \item{Parhist}{Matrix with paramater estimates across EM iterations}
#' \item{Mean}{M*(N*2) Matrix with NB means for every replicate and HMM component. The first two columns of Mean are the background and enrichment means of replicate 1, respectively, and so on}
#' \item{Viterbi}{Predicted sequence of Viterbi states}
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/ZIMHMM}
#'
#' @examples
#' data(H3K36me3.Huvec)
#' ChIP = as.matrix(H3K36me3.Huvec[,c("H3K36me3.Huvec.Rep1","H3K36me3.Huvec.Rep2","H3K36me3.Huvec.Rep3")])
#' Control = log(as.matrix(H3K36me3.Huvec[,c("Control.Huvec.Rep1","Control.Huvec.Rep2","Control.Huvec.Rep3")])+1)
#' offset = matrix(0,nrow = nrow(ChIP),ncol = ncol(ChIP),byrow = TRUE)
#' \dontrun{HMM(ChIP.init = rowSums(ChIP),Control.init = rowMeans(Control),offset.init = rowMeans(offset))}
#'
#' @export
#'
HMM = function(ChIP.init,Control.init,offset.init,pcut=0.05,epsilon.em=1e-3,maxit.em=5,minit.em=3,gap.em=3,maxcount.em=3,max.phi=1e3,min.zero=.Machine$double.xmin,quant=0.75,quiet=F){
    # Making data.table copies
    ChIP.init <- data.table::copy(ChIP.init)
    Control.init <- data.table::copy(Control.init)
    offset.init <- data.table::copy(offset.init)

    # General parameters
    M=length(ChIP.init);N=1;K=2
    error.em=1
    it.em = 0
    count.em = 0
    parlist = list()

    # Transforming data into data.table and calculating scores
    if(is.null(Control.init)){
        dt <- data.table(ChIP = ChIP.init,Dsg.Int = 1,offset = offset.init,PostProb1 = 1,PostProb2 = 1,JoinProb11 = 1,JoinProb12 = 1,JoinProb21 = 1,JoinProb22 = 1)
        Score = dt[,scale(log(ChIP+1))]
        ncolControl = 1
        namesControl = gsub('Dsg.','',names(dt)[grepl('Dsg',names(dt))])

        # Creating Aggragating Variable
        dt[,Group := .GRP,by=c('ChIP','Dsg.Int','offset')]

        # Creating Unique data.table
        dt.unique <- unique(dt,by='Group')[,c('ChIP','Dsg.Int','offset','Group')]
        setkey(dt.unique,Group)
    } else{
        dt <- data.table(ChIP = ChIP.init,Dsg.Int = 1,Dsg.Control = Control.init,offset = offset.init,PostProb1 = 1,PostProb2 = 1,JoinProb11 = 1,JoinProb12 = 1,JoinProb21 = 1,JoinProb22 = 1)
        Score = dt[,scale(log(ChIP+1)-Dsg.Control)]
        ncolControl = 2
        namesControl = gsub('Dsg.','',names(dt)[grepl('Dsg',names(dt))])

        # Creating Aggragating Variable
        dt[,Group := .GRP,by=c('ChIP','Dsg.Int','Dsg.Control','offset')]

        # Creating Unique data.table
        dt.unique <- unique(dt,by='Group')[,c('ChIP','Dsg.Int','Dsg.Control','offset','Group')]
        setkey(dt.unique,Group)
    }

    # Parameter initializations
    ## Hiden States
    dt[,z := as.numeric(cut(Score,breaks=c(-Inf,quantile(Score,quant),Inf)))]

    ##Initial probabilities
    pi1.old = 0.99;pi2.old = 1 - pi1.old;pi.old = c(pi1.old,pi2.old)

    ## Slope and Intercept of each component
    ### Aggregating data
    dt1 <- agg(data = dt,data.unique = dt.unique,rows = '(z==1)',agg = 'PostProb1')
    dt2 <- agg(data = dt,data.unique = dt.unique,rows = '(z==2)',agg = 'PostProb2')

    ### Calculating MLEs
    tryCatch({assign('model1i',optim(par=c(rep(0.25,ncolControl),1),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                     Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=F]),offset.vec=dt1[,offset],weights.vec=dt1[,weights]))},
             error=function(e){model1i<<-list();model1i[['par']]<<-c(rep(0.25,ncolControl),1)})
    tryCatch({assign('model2i',optim(par=c(rep(1,ncolControl),1),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                     Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=F]),offset.vec=dt2[,offset],weights.vec=dt2[,weights]))},
             error=function(e){model2i<<-list();model2i[['par']]<<-c(rep(1,ncolControl),1)})

    rm(dt1);rm(dt2)

    ### Saving parameters
    psi1.old = inv.par(model1i$par,model='nb');psi2.old = inv.par(model2i$par,model='nb');psi.old = c(psi1.old,psi2.old)

    ## Transition probabilities
    gamma.old = HMM.chain(dt[,z],K)

    # Putting all together
    theta.old = c(pi.old,gamma.old,psi.old)
    theta.k = theta.old
    names(theta.k) = c(paste0('pi',1:K),paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx)),
                       paste0('HMM1.',c(namesControl,'Disp')),paste0('HMM2.',c(namesControl,'Disp')))

    # EM algorithm begins
    while(count.em<maxcount.em & it.em<=maxit.em){
        it.em = it.em + 1

        # Updating parameters
        pi.k = theta.k[paste0('pi',1:K)]
        gamma.k = matrix(theta.k[paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx))],nrow=K,ncol=K,byrow=F);k=(K+1);for(i in 1:K){for(j in 1:K){assign(paste0('gamma',j,i,'.k'),theta.k[k]);k=k+1}}
        psi1.k = theta.k[paste0('HMM1.',c(namesControl,'Disp'))]
        psi2.k = theta.k[paste0('HMM2.',c(namesControl,'Disp'))]

        # E-step
        mu <- HMM.mean(X.mat=as.matrix(dt[,grepl('Dsg',names(dt)),with=F]),offset.vec=dt[,offset],psi=rbind(psi1.k,psi2.k)[,1:ncolControl],N=N,M=M)
        loglik <- HMM.LL(Y.vec=dt[,ChIP],mu=mu,disp=rbind(psi1.k,psi2.k)[,(ncolControl+1)],N=N,M=M,K=K,model='nb')

        # Forward-Backward probabilities
        logF <- hmm_logF(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)
        logB <- hmm_logB(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)

        # Posterior probabilities
        dt[,paste0('PostProb',1:2):=as.data.table(check.prob(hmm_P1(logF=logF,logB=logB)))]
        dt[,paste0('JoinProb',c('11','12','21','22')):=as.data.table(check.prob(hmm_P2(logF=logF,logB=logB,logf1=loglik[,1],logf2=loglik[,2],gamma=gamma.k)))]

        # M-step
        ## Initial and transition probabilities
        PostProb = HMM.prob(dt = dt)
        pi.k1 = PostProb$pi
        gamma.k1 = PostProb$gamma

        ## Model parameters
        ### Aggregating data
        rejection = (pcut>0)*ifelse((0.9^it.em)>=pcut,(0.9^it.em),pcut)
        dt1 <- agg(dt[,Rejection1 := PostProb1][PostProb1<rejection,Rejection1 := rbinom(.N,1,prob=PostProb1/rejection)*rejection],data.unique = dt.unique,rows = '(Rejection1>0)',agg = 'Rejection1')
        dt2 <- agg(dt[,Rejection2 := PostProb2][PostProb2<rejection,Rejection2 := rbinom(.N,1,prob=PostProb2/rejection)*rejection],data.unique = dt.unique,rows = '(Rejection2>0)',agg = 'Rejection2')

        ### Calculating MLEs
        tryCatch({assign('model1',optim(par=inv.par(psi1.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                        Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=F]),offset.vec=dt1[,offset],weights.vec=dt1[,weights]))},
                 error=function(e){model1<<-list();model1[['par']]<<-inv.par(psi1.k,model='nb');model1[['convergence']]<<-99})
        tryCatch({assign('model2',optim(par=inv.par(psi2.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                        Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=F]),offset.vec=dt2[,offset],weights.vec=dt2[,weights]))},
                 error=function(e){model2<<-list();model2[['par']]<<-inv.par(psi2.k,model='nb');model2[['convergence']]<<-99})

        rm(dt1);rm(dt2)

        ### Saving parameters
        psi1.k1 = inv.par(model1$par,model='nb')
        psi2.k1 = inv.par(model2$par,model='nb')
        psi.k1 = c(psi1.k1,psi2.k1)

        # Updating parameter history
        theta.k1 = c(pi.k1,gamma.k1,psi.k1)
        names(theta.k1) = names(theta.k)
        theta.k = theta.k1
        parlist[[it.em]] = c(it=it.em,error=error.em,theta.k1,m1=model1$convergence,m2=model2$convergence)

        # Computing EM error
        gap = ifelse(it.em>minit.em,gap.em,1)
        if(it.em>1){
            parlist.old = parlist[[(it.em-gap)]][names(psi.k1)]
            parlist.new = parlist[[it.em]][names(psi.k1)]
        } else{
            parlist.old = rep(1,length(names(psi.k1)))
            parlist.new = rep(1,length(names(psi.k1)))
        }
        MRCPE = max(abs((parlist.new-parlist.old)/parlist.old)) #Max. Abs. Rel. Change. of par. estimates
        error.em = ifelse(it.em>=2,MRCPE,1)
        count.em = as.numeric(any(error.em<=epsilon.em))*(it.em>minit.em)*(count.em+1) + 0
    }

    # Organizing output
    z = Viterbi(LOGF=loglik,P=pi.k1,GAMMA=gamma.k1)
    logF <- setnames(as.data.table(logF),c('Background','Enrichment'))
    logB <- setnames(as.data.table(logB),c('Background','Enrichment'))
    loglik <- setnames(as.data.table(loglik),c('Background','Enrichment'))
    mu <- as.data.table(mu)
    return(list('Pi'=pi.k1,'Gamma'=gamma.k1,'Psi'=psi.k1,'Prob'=dt[,.(PostProb1,PostProb2)],
                'LogF'=logF,'LogB'=logB,'Loglik'=loglik,'Parhist'=as.data.table(do.call(rbind,parlist)),'Mean'=mu,'Viterbi'=z))
}

