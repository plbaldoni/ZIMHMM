#' Zero Inflated Mixed Effects Hidden Markov Model (ZIMHMM)
#'
#' This function runs ZIMHMM.
#'
#' @param ChIP M*N matrix of ChIP read counts, where M is the number of windows in the analyzed genome and N is the number of replicates
#' @param Control M*N matrix of log-transformed Control read counts
#' @param offset M*N matrix of offsets. If no offset is used, use offset = matrix(0,nrow=M,ncol=N)
#' @param random either 'intercept' for random intercept model or 'slope' for random slope model
#' @param control list of control arguments from controlPeaks()
#'
#' @return A list with components
#' \item{Pi}{Vector of initial probabilities of the HMM}
#' \item{Gamma}{Matrix of transition probabilities of the HMM}
#' \item{Psi}{Vector of component-specific parameters of the HMM}
#' \item{Sigma2}{Variance component}
#' \item{U}{Vector of random effects}
#' \item{Zeroinfl}{M*N Matrix with zero-inflation probabilities}
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
#' data(Huvec)
#' ChIP = SummarizedExperiment::assay(Huvec,'ChIP')
#' Control = log(SummarizedExperiment::assay(Huvec,'Control')+1)
#' offset = matrix(0,nrow = nrow(ChIP),ncol = ncol(ChIP))
#' \dontrun{ZIMHMM(ChIP = ChIP,Control = Control,offset = offset,random = 'intercept',control = controlPeaks())}
#'
#' @export
#'
ZIMHMM <- function(ChIP,Control,offset,random,control){
    # Creating control elements
    for(i in seq_along(control)){assign(names(control)[i],control[[i]])}
    if(!(length(epsilon.em)==4) & criterion=='MULTI'){stop('For MULTI criterion, the length of error.em must be 4.')}

    # General parameters
    M=nrow(ChIP);N=ncol(ChIP);K=2
    error.em=1
    it.em = 0
    count.em = 0
    parlist = list()
    zlist = list()

    # Transforming data into data.table
    if(is.null(Control)){
        DT = data.table(ChIP = ChIP,Dsg.Int = 1,offset = offset,PostProb1=1,PostProb2=1,JoinProb11=1,JoinProb12=1,JoinProb21=1,JoinProb22=1,Rejection1=1,Rejection2=1)
        setnames(DT,c(paste0('ChIP.',1:N),'Dsg.Int',paste0('offset.',1:N),'PostProb1','PostProb2','JoinProb11','JoinProb12','JoinProb21','JoinProb22','Rejection1','Rejection2'))
        ncolControl = 1
        namesControl = c('Int')

        # Creating covariate associated with the random component (sigma2*u*1 if intercept, sigma2*u*Control if slope)
        ## Variance component and Latent vector of random effects
        sigma2.old = 0.10
        u.old = as.numeric(scale(colSums(log(ChIP+1))))

        # Creating covariate associates with the random component (1 if intercept, Control if slope)
        if(random=='intercept'){
            DT[,paste0('Random.',1:N) := as.data.table(DT[,mapply("*",sqrt(sigma2.old)*u.old,mget(rep('Dsg.Int',N)))])]
            DT[,paste0('U.',1:N) := as.list(u.old)]
            DT[,paste0('W.',1:N) := mget(rep('Dsg.Int',N))]
        } else{
            stop('Specify random="intercept"')
        }

        # Stacking DT
        DTvec <- vecData(DT,N,control = F,random = T)

        # Creating Aggragating Variable
        DTvec[,Group := .GRP,by=c('ChIP','Dsg.Int','offset','Random','U','W')]

        # Creating Unique data.table
        DTvec.unique <- unique(DTvec,by='Group')[,c('ChIP','Dsg.Int','offset','Random','U','W','Group')]
        setkey(DTvec.unique,Group)
    } else{
        DT = data.table(ChIP = ChIP,Dsg.Int = 1,Dsg.Control = Control,offset = offset,PostProb1=1,PostProb2=1,JoinProb11=1,JoinProb12=1,JoinProb21=1,JoinProb22=1,Rejection1=1,Rejection2=1)
        setnames(DT,c(paste0('ChIP.',1:N),'Dsg.Int',paste0('Dsg.Control.',1:N),paste0('offset.',1:N),'PostProb1','PostProb2','JoinProb11','JoinProb12','JoinProb21','JoinProb22','Rejection1','Rejection2'))
        ncolControl = 2
        namesControl = c('Int','Control')

        # Creating covariate associated with the random component (sigma2*u*1 if intercept, sigma2*u*Control if slope)
        ## Variance component and Latent vector of random effects
        sigma2.old = 0.10
        u.old = as.numeric(scale(colSums(log(ChIP+1))))

        if(random=='intercept'){
            DT[,paste0('Random.',1:N) := as.data.table(DT[,mapply("*",sqrt(sigma2.old)*u.old,mget(rep('Dsg.Int',N)))])]
            DT[,paste0('U.',1:N) := as.list(u.old)]
            DT[,paste0('W.',1:N) := mget(rep('Dsg.Int',N))]
        } else{
            DT[,paste0('Random.',1:N) := as.data.table(DT[,mapply("*",sqrt(sigma2.old)*u.old,mget(paste0('Dsg.Control.',1:N)))])]
            DT[,paste0('U.',1:N) := as.list(u.old)]
            DT[,paste0('W.',1:N) := mget(paste0('Dsg.Control.',1:N))]
        }

        # Stacking DT
        DTvec <- vecData(DT,N,random = T)

        # Creating Aggragating Variable
        DTvec[,Group := .GRP,by=c('ChIP','Dsg.Int','Dsg.Control','offset','Random','U','W')]

        # Creating Unique data.table
        DTvec.unique <- unique(DTvec,by='Group')[,c('ChIP','Dsg.Int','Dsg.Control','offset','Random','U','W','Group')]
        setkey(DTvec.unique,Group)
    }

    # Parameter initializations
    if(!quiet){cat(paste0(c(rep('#',45),'\n')));cat("Algorithm initialization...\n")}
    if(is.null(Control)){
        model.list = HMM(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
                              Control.init=NULL,
                              offset.init=DT[,rowMeans(.SD),.SDcols = paste0('offset.',1:N)],pcut=pcut)
    } else{
        model.list = HMM(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
                              Control.init=DT[,rowMeans(.SD),.SDcols = paste0('Dsg.Control.',1:N)],
                              offset.init=DT[,rowMeans(.SD),.SDcols = paste0('offset.',1:N)],pcut=pcut)
    }
    DT[,z := model.list$Viterbi]
    DTvec[,z := rep(model.list$Viterbi,N)]

    ## Initial probabilities
    pi1.old = 0.99;pi2.old = 1 - pi1.old;pi.old = c(pi1.old,pi2.old)

    ## Model-specific parameters
    ### Aggregating data
    dt1 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(z==0)',agg = 'PostProb1')
    dt2 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(z==1)',agg = 'PostProb2')

    ### Calculating MLEs
    tryCatch({assign('psi1.old',inv.par(optim(par=c(rep(0.25,ncolControl),1),fn=glm.nb,gr=deriv.nb,
                                              Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=F]),offset.vec=dt1[,offset],weights.vec=dt1[,weights],method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi))$par,model='nb')[c(1:ncolControl,(ncolControl+1),1:ncolControl)])},
             error=function(e){psi1.old<<-inv.par(c(rep(0.25,ncolControl),1),model='nb')[c(1:ncolControl,(ncolControl+1),1:ncolControl)]})
    tryCatch({assign('psi2.old',inv.par(optim(par=c(rep(1,ncolControl),1),fn=glm.nb,gr=deriv.nb,
                                              Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=F]),offset.vec=dt2[,offset],weights.vec=dt2[,weights],method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi))$par,model='nb'))},
             error=function(e){psi2.old<<-inv.par(c(rep(1,ncolControl),1),model='nb')})

    rm(dt1);rm(dt2)

    ### Saving parameters
    psi.old = c(psi1.old,psi2.old)

    ## Transition probabilities
    gamma.old = HMM.chain(DT[,z],K)

    #Putting all together
    theta.old = c(pi.old,gamma.old,psi.old,sigma2.old,u.old)
    theta.k = theta.old
    name.psi1 = c(paste0('HMM1.',c(namesControl,'Disp')),paste0('HMM1.ZI.',namesControl))
    names(theta.k) = c(paste0('pi',1:K),paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx)),
                       name.psi1,paste0('HMM2.',c(namesControl,'Disp')),'sigma2',paste0('U',1:N))

    if(!quiet){cat(paste0("Initialization completed!\n"));cat(paste0(c(rep('#',45),'\n')))}

    #EM algorithm begins
    if(!quiet){cat("EM algorithm begins...\n")}

    while(count.em<maxcount.em & it.em<=maxit.em){
        it.em = it.em+1

        #Updating parameters
        pi.k = theta.k[paste0('pi',1:K)]
        gamma.k = matrix(theta.k[paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx))],nrow=K,ncol=K,byrow=F);k=(K+1);for(i in 1:K){for(j in 1:K){assign(paste0('gamma',j,i,'.k'),theta.k[k]);k=k+1}}
        psi1.k = theta.k[c(paste0('HMM1.',c(namesControl,'Disp')),paste0('HMM1.ZI.',namesControl))]
        psi2.k = theta.k[paste0('HMM2.',c(namesControl,'Disp'))]
        sigma2.k = theta.k['sigma2']
        u.k = theta.k[paste0('U',1:N)]

        #E-step
        zeroinfl <- HMM.zeroinfl(csi=psi1.k[(ncolControl+2):length(psi1.k)],X.mat=as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),offset.vec=DTvec[,offset],N=N,M=M)

        ## Laplace Approximation: 'estimating' latent random effects
        U.opt <- minqa::bobyqa(par = u.k,fn = function(x,...){-integrand(x,...)},
                               YVEC=DTvec[,ChIP],XMAT=as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),
                               BETA=matrix(c(psi1.k[1:ncolControl],psi2.k[1:ncolControl]),nrow=K,byrow=T),
                               DISP=c(psi1.k[(ncolControl+1)],psi2.k[(ncolControl+1)]),P=pi.k,GAMMA=gamma.k,
                               OFFSETVEC=DTvec[,offset],ZEROINFL=zeroinfl,
                               W=DTvec[,W],SIGMA2=sigma2.k)
        u.k1 = U.opt$par

        ########################################################################
        ## Updating DT and DTvec
        if(random=='intercept'){
            DT[,paste0('Random.',1:N) := as.data.table(DT[,mapply("*",sqrt(sigma2.k)*u.k1,mget(rep('Dsg.Int',N)))])]
            DT[,paste0('U.',1:N) := as.list(u.k1)]
            DT[,paste0('W.',1:N) := mget(rep('Dsg.Int',N))]
        } else{
            DT[,paste0('Random.',1:N) := as.data.table(DT[,mapply("*",sqrt(sigma2.k)*u.k1,mget(paste0('Dsg.Control.',1:N)))])]
            DT[,paste0('U.',1:N) := as.list(u.k1)]
            DT[,paste0('W.',1:N) := mget(paste0('Dsg.Control.',1:N))]
        }

        # Stacking DT
        DTvec <- vecData(DT,N,random = T)

        # Creating Aggragating Variable
        DTvec[,Group := .GRP,by=c('ChIP','Dsg.Int','Dsg.Control','offset','Random','U','W')]

        # Creating Unique data.table
        DTvec.unique <- unique(DTvec,by='Group')[,c('ChIP','Dsg.Int','Dsg.Control','offset','Random','U','W','Group')]
        setkey(DTvec.unique,Group)
        ########################################################################

        ## Mean and log-likelihood matrices
        mu = MHMMmean(XMAT = as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),BETA = matrix(c(psi1.k[1:ncolControl],psi2.k[1:ncolControl]),nrow=K,byrow=T),RANDOM = as.matrix(DTvec[,Random]),OFFSETVEC = as.matrix(DTvec[,offset]),N = N,M = M,K = K)
        mu[mu==0] = min.zero
        loglik = MHMMLik(YVEC = DTvec[,ChIP],ZEROINFL = zeroinfl,MU = mu,DISP = c(psi1.k[(ncolControl+1)],psi2.k[(ncolControl+1)]),N = N,M = M,K = K)

        # Forward-Backward probabilities
        logF = hmm_logF(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)
        logB = hmm_logB(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)

        # Posterior probabilities
        DT[,paste0('PostProb',1:2):=as.data.table(check.prob(hmm_P1(logF=logF,logB=logB)))]
        DT[,paste0('JoinProb',c('11','12','21','22')):=as.data.table(check.prob(hmm_P2(logF=logF,logB=logB,logf1=loglik[,1],logf2=loglik[,2],gamma=gamma.k)))]

        # M-step
        ## Initial and transition probabilities
        PostProb = HMM.prob(DT)
        pi.k1 = PostProb$pi
        gamma.k1 = PostProb$gamma
        zlist[[it.em]] = Viterbi(LOGF=loglik,P=pi.k1,GAMMA=gamma.k1)

        ## Model parameters
        ### General parameters for conditional EM iterations
        error.inner.em = 1
        count.inner.em = 1
        parhist.inner.em = data.frame()

        ### Updating posterior probabilities with rejection-controlled EM
        rejection = (pcut>0)*ifelse((0.9^it.em)>=pcut,(0.9^it.em),pcut)
        DT[,c('Rejection1','Rejection2') := list(PostProb1,PostProb2)]
        DT[PostProb1<rejection,Rejection1 := rbinom(.N,1,prob=PostProb1/rejection)*rejection]
        DT[PostProb2<rejection,Rejection2 := rbinom(.N,1,prob=PostProb2/rejection)*rejection]

        ### Updating the vectorized dataset
        DTvec[,c('Rejection1','Rejection2') := .(rep(DT[,Rejection1],N),rep(DT[,Rejection2],N))]
        DTvec[,c('PostProb1','PostProb2') := .(rep(DT[,PostProb1],N),rep(DT[,PostProb2],N))]

        ### Aggregating data
        dt1 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(Rejection1>0)',agg = 'Rejection1')
        dt2 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(Rejection2>0)',agg = 'Rejection2')
        dtsigma <- agg(data = DTvec,data.unique = DTvec.unique,agg = c('PostProb1','PostProb2'),random = T)

        ### Conditional EM begins
        while(error.inner.em>epsilon.inner.em & count.inner.em<=maxcount.inner.em){
            tryCatch({assign('model1',optim(par = inv.par(psi1.k,'zinb'),fn = glmm.zinb,gr = deriv.glmm.zinb,method = 'L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi,rep(-Inf,ncolControl)),
                                            Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=F]),offset.glm.vec=as.matrix(dt1[,rowSums(.SD),.SDcols=c('Random','offset')]),offset.zi.vec=as.matrix(dt1[,offset]),weights.vec=as.matrix(dt1[,weights])))},
                     error=function(e){model1<<-list();model1[['par']]<<-inv.par(psi1.k,'zinb');model1[['convergence']]<<-99})
            tryCatch({assign('model2',optim(par=inv.par(psi2.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                            Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=F]),offset.vec=as.matrix(dt2[,rowSums(.SD),.SDcols=c('Random','offset')]),weights.vec=as.matrix(dt2[,weights])))},
                     error=function(e){model2<<-list();model2[['par']]<<-inv.par(psi2.k,model='nb');model2[['convergence']]<<-99})

            psi1.k = inv.par(model1$par,model='zinb');names(psi1.k) = name.psi1
            psi2.k = inv.par(model2$par,model='nb')

            tryCatch({assign('model3',optim(par = 1/sigma2.old,fn = glm.s2,method = 'L-BFGS-B',lower = 1/max.sigma2,upper = 1/min.sigma2,
                                            Yvec = dtsigma[,ChIP],Xmat = as.matrix(dtsigma[,grepl('Dsg',names(dtsigma)),with=F]),Uvec = dtsigma[,U],Wvec = dtsigma[,W],
                                            PostProbBackground = dtsigma[,AggPostProb1],PostProbEnrichment = dtsigma[,AggPostProb2],Offset = dtsigma[,offset],
                                            BackgroundPar = psi1.k,EnrichmentPar = psi2.k))},
                     error=function(e){model3<<-list();model3[['par']]<<-1/sigma2.k;model3[['convergence']]<<-99})

            sigma2.k = 1/model3$par

            # Updating the aggregated data
            dt1[,Random := sqrt(sigma2.k)*U]
            dt2[,Random := sqrt(sigma2.k)*U]

            # Updating inner parameters
            parhist.inner.em = rbind(parhist.inner.em,data.frame(t(c(psi1.k,psi2.k,sigma2=sigma2.k))))

            # Calculating inner error
            if(count.inner.em>1){error.inner.em = max(abs((parhist.inner.em[count.inner.em,]-parhist.inner.em[(count.inner.em-1),])/parhist.inner.em[(count.inner.em-1),]))}
            count.inner.em = count.inner.em + 1
        }

        # Updating parameter history
        psi1.k1 = psi1.k;psi2.k1 = psi2.k;psi.k1 = c(psi1.k1,psi2.k1)
        sigma2.k1 = sigma2.k
        theta.k1 = c(pi.k1,gamma.k1,psi.k1,sigma2.k1,u.k1);names(theta.k1) = names(theta.k)
        theta.k = theta.k1
        parlist[[it.em]] = c(it=it.em,Q=Q(as.matrix(DT[,.(PostProb1,PostProb2)]),as.matrix(DT[,.(JoinProb11,JoinProb12,JoinProb21,JoinProb22)]),loglik,pi.k1,gamma.k1),
                             error=error.em[1],theta.k1,m1=model1$convergence,m2=model2$convergence,m3=model3$convergence)

        # Computing EM error
        gap = ifelse(it.em>minit.em,gap.em,1)
        if(it.em>1){
            parlist.old = parlist[[(it.em-gap)]][names(psi.k1)]
            parlist.new = parlist[[it.em]][names(psi.k1)]
            zlist.table = data.table(old = zlist[[(it.em-gap)]], new = zlist[[it.em]])
            ACC = 100*zlist.table[,.N,by=.(old,new)][(old==0 & new==0) | (old==1 & new==1),sum(N)]/M
        } else{
            parlist.old = rep(1,length(names(psi.k1)))
            parlist.new = rep(1,length(names(psi.k1)))
            ACC = 0
        }

        MRCPE = max(abs((parlist.new-parlist.old)/parlist.old)) #Max. Abs. Rel. Change. of par. estimates
        MACPE = max(abs(parlist.new-parlist.old)) #Max. Abs. Change. of par. estimates
        ARCEL = ifelse(it.em>=2,abs((parlist[[it.em]][['Q']] - parlist[[(it.em-gap)]][['Q']])/parlist[[(it.em-gap)]][['Q']]),0) #Abs. Rel. Change of expected log-likelihood of complete data (Q function)
        MULTI = c(MRCPE,MACPE,ARCEL,100-ACC)
        error.em = (it.em>=2)*get(criterion) + (it.em<2)*rep(1,length(get(criterion)))
        count.em = as.numeric(any(error.em<=epsilon.em))*(it.em>minit.em)*(count.em+1) + 0

        #Outputing history
        if(!quiet){
            cat(paste0(c(rep('#',45),'\n')))
            cat('\rIteration: ',it.em,', Error(s): ',paste(formatC(error.em, format = "e", digits = 2),collapse = ', '),', Viterbi Agreement: ',round(ACC,2),'%.\n',sep='')
            cat("\r",paste('Q-function: '),parlist[[it.em]][['Q']],"\n")
            cat("\r",paste('Max. abs. rel. change of parameter estimates: '),MRCPE,"\n")
            cat("\r",paste('Max. abs. change of parameter estimates: '),MACPE,"\n")
            cat("\r",paste('Abs. rel. change of Q-function: '),ARCEL,"\n")
            cat(paste0(c(rep('#',45),'\n')))
        }
    }

    # Organizing output
    logF <- setnames(as.data.table(logF),c('Background','Enrichment'))
    logB <- setnames(as.data.table(logB),c('Background','Enrichment'))
    loglik <- setnames(as.data.table(loglik),c('Background','Enrichment'))
    mu <- as.data.table(mu)
    zeroinfl <- as.data.table(HMM.zeroinfl(csi=psi1.k[(ncolControl+2):length(psi1.k)],X.mat=as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),offset.vec=DTvec[,offset],N=N,M=M))

    if(!quiet){cat('\nDone!\n')}
    return(list('Pi'=pi.k1,'Gamma'=gamma.k1,'Psi'=psi.k1,'Sigma2'=sigma2.k1,'U'=u.k1,
                'Zeroinfl'=zeroinfl,'Prob'=DT[,.(PostProb1,PostProb2)],'LogF'=logF,'LogB'=logB,'Loglik'=loglik,
                'Parhist'=as.data.frame(do.call(rbind,parlist)),'Mean'=mu,'Viterbi'=zlist[[it.em]]))
}
