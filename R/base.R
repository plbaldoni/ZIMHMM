vecData = function(data,N,control=T,random=F){
    data <- data.table::copy(data)
    if(!random){
        if(control){
            return(melt(data[,unique(c(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('Dsg.Control.',1:sum(grepl('^Dsg.Control',names(data)))),paste0('offset.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N))),with=F],
                        measure = list(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('Dsg.Control.',1:sum(grepl('^Dsg.Control',names(data)))),paste0('offset.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N)),
                        value.name = c('ChIP','Dsg.Int','Dsg.Control','offset','PostProb1','PostProb2','Rejection1','Rejection2'))[,-1])
        } else{
            return(melt(data[,unique(c(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('offset.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N))),with=F],
                        measure = list(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('offset.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N)),
                        value.name = c('ChIP','Dsg.Int','offset','PostProb1','PostProb2','Rejection1','Rejection2'))[,-1])
        }
    } else{
        if(control){
            return(melt(data[,unique(c(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('Dsg.Control.',1:sum(grepl('^Dsg.Control',names(data)))),paste0('offset.',1:N),paste0('Random.',1:N),paste0('U.',1:N),paste0('W.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N))),with=F],
                        measure = list(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('Dsg.Control.',1:sum(grepl('^Dsg.Control',names(data)))),paste0('offset.',1:N),paste0('Random.',1:N),paste0('U.',1:N),paste0('W.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N)),
                        value.name = c('ChIP','Dsg.Int','Dsg.Control','offset','Random','U','W','PostProb1','PostProb2','Rejection1','Rejection2'))[,-1])
        } else{
            return(melt(data[,unique(c(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('offset.',1:N),paste0('Random.',1:N),paste0('U.',1:N),paste0('W.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N))),with=F],
                        measure = list(paste0('ChIP.',1:N),rep('Dsg.Int',N),paste0('offset.',1:N),paste0('Random.',1:N),paste0('U.',1:N),paste0('W.',1:N),rep('PostProb1',N),rep('PostProb2',N),rep('Rejection1',N),rep('Rejection2',N)),
                        value.name = c('ChIP','Dsg.Int','offset','Random','U','W','PostProb1','PostProb2','Rejection1','Rejection2'))[,-1])
        }
    }
}

agg = function(data,data.unique,rows,agg,random=F){
    data <- data.table::copy(data)
    setkey(data,Group)

    data.unique <- data.table::copy(data.unique)

    # Filtering rows and aggregating
    # This will give me the sum of weights per each Group (possible combination of ChIP, Control, etc.)
    if(random==F){
        data <- data[eval(parse(text=rows)),][,.(weights=sum(get(agg))),by=Group]
    }
    if(random==T){
        data <- data[,.(AggPostProb1=sum(get(agg[1])),AggPostProb2=sum(get(agg[2]))),by=Group]
    }

    # Bringing the variables from data.unique into the aggregated data
    # This will bring the variables back to the data
    return(data.unique[data,on='Group'])

    # return(data[eval(parse(text=rows)),][,.(weights=sum(get(agg))),by=eval(names(data)[!grepl('^PostProb|^JoinProb|^Rejection',names(data))])])
}

createOffset = function(ChIP,method='sum',span=0.3,plots=F,dirplot=NULL){
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

Q = function(P1,P2,loglik,pi,gamma){
    return(sum(P1[1,]*log(pi))+sum(colSums(P1*loglik))+sum(colSums(P2[2:nrow(loglik),]%*%diag(log(c(t(gamma)))))))
}

glm.nb = function(par,Y.vec,X.mat,offset.vec,weights.vec){
    l = dnbinom(Y.vec,mu=exp(X.mat%*%par[1:ncol(X.mat)]+offset.vec),size=1/par[length(par)],log=T);l[is.infinite(l)] = log(1e-300)
    return(-sum(weights.vec*l))
}

deriv.nb = function(par,Y.vec,X.mat,offset.vec,weights.vec,...){
    mu.vec = exp(X.mat%*%par[1:ncol(X.mat)]+offset.vec)
    phi = par[length(par)]
    return(-c(colSums(as.numeric(weights.vec)*as.numeric((Y.vec-mu.vec)/(1+phi*mu.vec))*X.mat),sum(as.numeric(weights.vec)*(log(1+phi*mu.vec)+phi*(Y.vec-mu.vec)/(1+phi*mu.vec)-digamma(Y.vec+1/phi)+digamma(1/phi))/(phi^2))))
}

inv.par = function(par,model){
    if(model=='nb'){return(c(par[1:(length(par)-1)],1/par[length(par)]))}
    if(model=='zinb'){
        n = length(par)
        return(c(par[1:((n-1)/2)],1/par[(n-1)/2+1],par[((n-1)/2+2):n]))}
}

HMM.chain = function(z,K){
    if(max(z)>(K-1)){z=round((K-1)*(z-max(z))/(max(z)-min(z))+(K-1))}
    MC = matrix(z,nrow = 1,ncol=length(z))
    MC = table(c(MC[,-ncol(MC)]),c(MC[,-1]))
    MC = as.matrix(MC/rowSums(MC))
    MC = matrix(MC,ncol=ncol(MC),nrow=nrow(MC),byrow=F)
    MC = check.prob(MC)
    return(MC)
}

HMM.mean = function(X.mat,offset.vec,psi,N,M,min.zero=min.zero,U=NULL,random=NULL){
    K=2#nrow(psi)
    mu = matrix(0,nrow=M,ncol=N*K)
    for(k in 1:K){
        mu[,k+0:(N-1)*K] = matrix(exp(X.mat%*%psi[k,]+offset.vec),nrow=M,ncol=N,byrow=F)
    }
    if(sum(mu==0)>0){mu[mu==0] = min.zero}
    return(mu)
}

HMM.LL = function(Y.vec,mu,N,M,K,model,disp=NULL,zeroinfl=NULL,min.zero=.Machine$double.xmin){
    LL = matrix(0,nrow=M,ncol=K)
    # if(model=='poisson'){
    #     for(k in 1:K){
    #         LL[,k] = rowSums(matrix((dpois(Y.vec,lambda=c(mu[,k+0:(N-1)*K]),log=T)),nrow=M,ncol=N,byrow=F))
    #     }
    #
    #     LL[exp(LL)==0] = log(min.zero)
    #     return(LL)
    # }
    if(model=='nb'){
        for(k in 1:K){
            LL[,k] = rowSums(matrix((dnbinom(Y.vec,mu=c(mu[,k+0:(N-1)*K]),size=disp[k],log=T)),nrow=M,ncol=N,byrow=F))
        }
        LL[is.infinite(LL)] = log(min.zero)
        # LL[exp(LL)==0] = log(min.zero)
        return(LL)
    }
    # if(model=='zip'){
    #     Yvec0 = 1*(Y.vec==0)
    #     Zvec = c(zeroinfl)
    #
    #     #LL from Background
    #     LogLik00 = dpois(0,lambda=c(mu[,1+0:(N-1)*K]),log=T);LogLik00[is.infinite(LogLik00)] = log(1e-300)
    #     LogLik01 = dpois(Y.vec,lambda=c(mu[,1+0:(N-1)*K]),log=T);LogLik01[is.infinite(LogLik01)] = log(1e-300)
    #     LL[,1] = rowSums(matrix((Yvec0)*log(Zvec + exp(log(1-Zvec)+LogLik00)) + (1-Yvec0)*(log(1-Zvec) + LogLik01),nrow=M,ncol=N,byrow=F))
    #     LL[is.infinite(LL[,1]),1] = log(1e-300)
    #     #LL from Differential
    #     LL[,2] = rowSums(matrix((dpois(Y.vec,lambda=c(mu[,2+0:(N-1)*K]),log=T)),nrow=M,ncol=N,byrow=F))
    #     #LL from Enrichment
    #     LL[,K] = rowSums(matrix((dpois(Y.vec,lambda=c(mu[,K+0:(N-1)*K]),log=T)),nrow=M,ncol=N,byrow=F))
    #
    #     LL[exp(LL)==0] = log(min.zero)
    #     return(LL)
    # }
    if(model=='zinb'){
        Yvec0 = 1*(Y.vec==0)
        Zvec = c(zeroinfl)
        # if(K==2){
            #LL from Background
            LogLik00 = dnbinom(0,mu=c(mu[,1+0:(N-1)*K]),size=disp[1],log=T);LogLik00[is.infinite(LogLik00)] = log(min.zero)
            LogLik01 = dnbinom(Y.vec,mu=c(mu[,1+0:(N-1)*K]),size=disp[1],log=T);LogLik01[is.infinite(LogLik01)] = log(min.zero)
            LL[,1] = rowSums(matrix((Yvec0)*log(Zvec + exp(log(1-Zvec)+LogLik00)) + (1-Yvec0)*(log(1-Zvec) + LogLik01),nrow=M,ncol=N,byrow=F))
            LL[is.infinite(LL[,1]),1] = log(min.zero)

            #LL from Enrichment
            LL[,K] = rowSums(matrix((dnbinom(Y.vec,mu=c(mu[,K+0:(N-1)*K]),size=disp[K],log=T)),nrow=M,ncol=N,byrow=F))
            LL[is.infinite(LL[,K]),K] = log(min.zero)

            # LL[exp(LL)==0] = log(min.zero)
            return(LL)
        # }
        # if(K==3){
        #     #LL from Background
        #     LogLik00 = dnbinom(0,mu=c(mu[,1+0:(N-1)*K]),size=disp[1],log=T);LogLik00[is.infinite(LogLik00)] = log(1e-300)
        #     LogLik01 = dnbinom(Y.vec,mu=c(mu[,1+0:(N-1)*K]),size=disp[1],log=T);LogLik01[is.infinite(LogLik01)] = log(1e-300)
        #     LL[,1] = rowSums(matrix((Yvec0)*log(Zvec + exp(log(1-Zvec)+LogLik00)) + (1-Yvec0)*(log(1-Zvec) + LogLik01),nrow=M,ncol=N,byrow=F))
        #     LL[is.infinite(LL[,1]),1] = log(1e-300)
        #     #LL from Differential
        #     LL[,2] = rowSums(matrix((dnbinom(Y.vec,mu=c(mu[,2+0:(N-1)*K]),size=disp[2],log=T)),nrow=M,ncol=N,byrow=F))
        #     #LL from Enrichment
        #     LL[,K] = rowSums(matrix((dnbinom(Y.vec,mu=c(mu[,K+0:(N-1)*K]),size=disp[K],log=T)),nrow=M,ncol=N,byrow=F))
        #
        #     LL[exp(LL)==0] = log(min.zero)
        #     return(LL)
        # }
    }
}

check.prob = function(P){
    P = pmax(pmin(P,1),0)
    return(P/rowSums(P))
}

HMM.prob = function(dt){
    pi <- dt[1,c(PostProb1,PostProb2)]
    gamma <- unname(check.prob(as.matrix(t(dt[,.(c(sum(JoinProb11,na.rm=T),sum(JoinProb12,na.rm=T))/sum(JoinProb11+JoinProb12,na.rm=T),
                                                c(sum(JoinProb21,na.rm=T),sum(JoinProb22,na.rm=T))/sum(JoinProb21+JoinProb22,na.rm=T))]))))
    return(list('pi'=pi,'gamma'=gamma))
}

glm.zinb = function(par,Y.vec,X.mat,offset.vec,weights.vec){
    mu.vec = exp(X.mat%*%par[1:ncol(X.mat)]+offset.vec)
    zeroinfl.vec = 1/(1+exp(-(X.mat%*%par[(ncol(X.mat)+2):(2*ncol(X.mat)+1)]+offset.vec)))
    idx.Y0=(Y.vec==0)

    l0 = suppressWarnings(dnbinom(0,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l0[is.infinite(l0)] = log(1e-300)
    l1 = suppressWarnings(dnbinom(Y.vec,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l1[is.infinite(l1)] = log(1e-300)

    Loglik = weights.vec*(idx.Y0*log(zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0))+(1-idx.Y0)*(log(1-zeroinfl.vec)+l1))

    return(-sum(Loglik))
}

deriv.glm.zinb = function(par,Y.vec,X.mat,offset.vec,weights.vec){
    # This function is correct (Checked with grad() from NumDeriv on 11/03/18)
    mu.vec = exp(X.mat%*%par[1:ncol(X.mat)]+offset.vec)
    zeroinfl.vec = 1/(1+exp(-(X.mat%*%par[(ncol(X.mat)+2):(2*ncol(X.mat)+1)]+offset.vec)))
    idx.Y0=(Y.vec==0)

    l0 = suppressWarnings(dnbinom(0,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l0[is.infinite(l0)] = log(1e-300)
    l1 = suppressWarnings(dnbinom(Y.vec,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l1[is.infinite(l1)] = log(1e-300)

    P0 = zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0)
    P1 = exp(log(1-zeroinfl.vec)+l1)

    phi = par[(ncol(X.mat)+1)]
    aux = (1+phi*mu.vec)

    return(c(-colSums(as.numeric((idx.Y0)*(weights.vec/P0)*(1-zeroinfl.vec)*(-mu.vec*((1+phi*mu.vec)^(-(1/phi+1)))))*X.mat+as.numeric((1-idx.Y0)*(weights.vec/P1)*(1-zeroinfl.vec)*exp(l1)*(Y.vec-mu.vec)/(1+phi*mu.vec))*X.mat),
             -sum((idx.Y0)*(weights.vec/P0)*(1-zeroinfl.vec)*(aux^(-1/phi))*(log(aux)/((phi)^2)-mu.vec/(phi*aux))+(1-idx.Y0)*(weights.vec/P1)*(1-zeroinfl.vec)*exp(l1)*(log(aux)/(phi^2)-mu.vec*(Y.vec+1/phi)/aux-digamma(Y.vec+1/phi)/(phi^2)+digamma(1/phi)/(phi^2)+Y.vec/phi)),
             -colSums(as.numeric((idx.Y0)*(weights.vec/P0)*(1-exp(l0))*zeroinfl.vec*(1-zeroinfl.vec))*X.mat+as.numeric((1-idx.Y0)*(weights.vec/P1)*(-exp(l1))*zeroinfl.vec*(1-zeroinfl.vec))*X.mat)))
}

glmm.zinb = function(par,Y.vec,X.mat,offset.glm.vec,offset.zi.vec,weights.vec){
    mu.vec = exp(X.mat%*%par[1:ncol(X.mat)]+offset.glm.vec)
    zeroinfl.vec = 1/(1+exp(-(X.mat%*%par[(ncol(X.mat)+2):(2*ncol(X.mat)+1)]+offset.zi.vec)))
    idx.Y0=(Y.vec==0)

    l0 = suppressWarnings(dnbinom(0,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l0[is.infinite(l0)] = log(1e-300)
    l1 = suppressWarnings(dnbinom(Y.vec,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l1[is.infinite(l1)] = log(1e-300)

    Loglik = weights.vec*(idx.Y0*log(zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0))+(1-idx.Y0)*(log(1-zeroinfl.vec)+l1))

    return(-sum(Loglik))
}

deriv.glmm.zinb = function(par,Y.vec,X.mat,offset.glm.vec,offset.zi.vec,weights.vec){
    mu.vec = exp(X.mat%*%par[1:ncol(X.mat)]+offset.glm.vec)
    zeroinfl.vec = 1/(1+exp(-(X.mat%*%par[(ncol(X.mat)+2):(2*ncol(X.mat)+1)]+offset.zi.vec)))
    idx.Y0=(Y.vec==0)

    l0 = suppressWarnings(dnbinom(0,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l0[is.infinite(l0)] = log(1e-300)
    l1 = suppressWarnings(dnbinom(Y.vec,mu=mu.vec,size=1/par[(ncol(X.mat)+1)],log=T));l1[is.infinite(l1)] = log(1e-300)

    P0 = zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0)
    P1 = exp(log(1-zeroinfl.vec)+l1)

    phi = par[(ncol(X.mat)+1)]
    aux = (1+phi*mu.vec)

    return(c(-colSums(as.numeric((idx.Y0)*(weights.vec/P0)*(1-zeroinfl.vec)*(-mu.vec*((1+phi*mu.vec)^(-(1/phi+1)))))*X.mat+as.numeric((1-idx.Y0)*(weights.vec/P1)*(1-zeroinfl.vec)*exp(l1)*(Y.vec-mu.vec)/(1+phi*mu.vec))*X.mat),
             -sum((idx.Y0)*(weights.vec/P0)*(1-zeroinfl.vec)*(aux^(-1/phi))*(log(aux)/((phi)^2)-mu.vec/(phi*aux))+(1-idx.Y0)*(weights.vec/P1)*(1-zeroinfl.vec)*exp(l1)*(log(aux)/(phi^2)-mu.vec*(Y.vec+1/phi)/aux-digamma(Y.vec+1/phi)/(phi^2)+digamma(1/phi)/(phi^2)+Y.vec/phi)),
             -colSums(as.numeric((idx.Y0)*(weights.vec/P0)*(1-exp(l0))*zeroinfl.vec*(1-zeroinfl.vec))*X.mat+as.numeric((1-idx.Y0)*(weights.vec/P1)*(-exp(l1))*zeroinfl.vec*(1-zeroinfl.vec))*X.mat)))
}

HMM.zeroinfl = function(csi,X.mat,offset.vec,N,M){
    Z = matrix(1/(1+exp(-(X.mat%*%csi+offset.vec))),nrow=M,ncol=N)
    return(Z)
}

glm.s2 = function(par,
                  Yvec,Xmat,Uvec,Wvec,PostProbBackground,PostProbEnrichment,Offset,
                  BackgroundPar,EnrichmentPar,min.zero=.Machine$double.xmin){
    ### Zero indicator to be used in the ZINB model
    Yvec0 = 1*(Yvec==0)
    ### ZI model
    zeroinfl.vec = 1/(1+exp(-(Xmat%*%BackgroundPar[(ncol(Xmat)+2):(2*ncol(Xmat)+1)]+Offset)))

    ### Mean vectors
    mu0 = exp(Xmat%*%BackgroundPar[1:ncol(Xmat)]+sqrt(1/par)*Uvec*Wvec+Offset)
    mu1 = exp(Xmat%*%EnrichmentPar[1:ncol(Xmat)]+sqrt(1/par)*Uvec*Wvec+Offset)

    ### Background log-likelihood
    f00 = suppressWarnings(dnbinom(0,mu=mu0,size=BackgroundPar[(ncol(Xmat)+1)],log=F));f00[f00==0] = min.zero
    f01 = suppressWarnings(dnbinom(Yvec,mu=mu0,size=BackgroundPar[(ncol(Xmat)+1)],log=F));f01[f01==0] = min.zero
    f0 = log(Yvec0*(zeroinfl.vec+(1-zeroinfl.vec)*f00) + (1-Yvec0)*((1-zeroinfl.vec)*f01));f0[is.infinite(f0)] = log(min.zero)

    ### Enrichment log-likelihood
    f1 = suppressWarnings(dnbinom(Yvec,mu=mu1,size=EnrichmentPar[(ncol(Xmat)+1)],log=T));f1[is.infinite(f1)] = log(min.zero)

    return(-sum(PostProbBackground*f0)-sum(PostProbEnrichment*f1))
}

HMM.init = function(ChIP.init,Control.init,offset.init,pcut,epsilon.em=1e-3,maxit.em=5,minit.em=3,gap.em=3,maxcount.em=3,max.phi=1e3,min.zero=.Machine$double.xmin,quant=0.75,quiet=T){
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

ZIHMM = function(ChIP,Control,offset,control)
{
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

        # Stacking DT
        DTvec <- vecData(DT,N,control = F)

        # Creating Aggragating Variable
        DTvec[,Group := .GRP,by=c('ChIP','Dsg.Int','offset')]

        # Creating Unique data.table
        DTvec.unique <- unique(DTvec,by='Group')[,c('ChIP','Dsg.Int','offset','Group')]
        setkey(DTvec.unique,Group)
    } else{
        DT = data.table(ChIP = ChIP,Dsg.Int = 1,Dsg.Control = Control,offset = offset,PostProb1=1,PostProb2=1,JoinProb11=1,JoinProb12=1,JoinProb21=1,JoinProb22=1,Rejection1=1,Rejection2=1)
        setnames(DT,c(paste0('ChIP.',1:N),'Dsg.Int',paste0('Dsg.Control.',1:N),paste0('offset.',1:N),'PostProb1','PostProb2','JoinProb11','JoinProb12','JoinProb21','JoinProb22','Rejection1','Rejection2'))
        ncolControl = 2
        namesControl = c('Int','Control')

        # Stacking DT
        DTvec <- vecData(DT,N)

        # Creating Aggragating Variable
        DTvec[,Group := .GRP,by=c('ChIP','Dsg.Int','Dsg.Control','offset')]

        # Creating Unique data.table
        DTvec.unique <- unique(DTvec,by='Group')[,c('ChIP','Dsg.Int','Dsg.Control','offset','Group')]
        setkey(DTvec.unique,Group)
    }

    # Parameter initializations
    if(!quiet){cat(paste0(c(rep('#',45),'\n')));cat("Algorithm initialization...\n")}
    if(is.null(Control)){
        model.list = HMM.init(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
                              Control=NULL,
                              offset.init=DT[,rowMeans(.SD),.SDcols = paste0('offset.',1:N)],pcut=pcut)
    } else{
        model.list = HMM.init(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
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

    # Putting all together
    theta.old = c(pi.old,gamma.old,psi.old)
    theta.k = theta.old
    names(theta.k) = c(paste0('pi',1:K),paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx)),
                       c(paste0('HMM1.',c(namesControl,'Disp')),paste0('HMM1.ZI.',namesControl)),paste0('HMM2.',c(namesControl,'Disp')))

    if(!quiet){cat(paste0("Initialization completed!\n"));cat(paste0(c(rep('#',45),'\n')))}

    # EM algorithm begins
    if(!quiet){cat("EM algorithm begins...\n")}

    while(count.em<maxcount.em & it.em<=maxit.em){
        it.em = it.em+1

        # Updating parameters
        pi.k = theta.k[paste0('pi',1:K)]
        gamma.k = matrix(theta.k[paste0('gamma',as.character(transform(expand.grid(1:K,1:K),idx=paste0(Var1,Var2))$idx))],nrow=K,ncol=K,byrow=F);k=(K+1);for(i in 1:K){for(j in 1:K){assign(paste0('gamma',j,i,'.k'),theta.k[k]);k=k+1}}
        psi1.k = theta.k[c(paste0('HMM1.',c(namesControl,'Disp')),paste0('HMM1.ZI.',namesControl))]
        psi2.k = theta.k[paste0('HMM2.',c(namesControl,'Disp'))]

        # E-step
        mu <- HMM.mean(X.mat=as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),offset.vec=DTvec[,offset],psi=rbind(psi1.k[1:ncolControl],psi2.k[1:ncolControl]),N=N,M=M)
        zeroinfl <- HMM.zeroinfl(csi=psi1.k[(ncolControl+2):length(psi1.k)],X.mat=as.matrix(DTvec[,grepl('Dsg',names(DTvec)),with=F]),offset.vec=DTvec[,offset],N=N,M=M)
        loglik <- HMM.LL(Y.vec=DTvec[,ChIP],mu=mu,N=N,M=M,K=K,zeroinfl=zeroinfl,disp=c(psi1.k[ncolControl+1],psi2.k[ncolControl+1]),model='zinb')

        # Forward-Backward probabilities
        logF <- hmm_logF(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)
        logB <- hmm_logB(logf1 = loglik[,1], logf2 = loglik[,2], pi = pi.k,gamma=gamma.k)

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
        ### Updating posterior probabilities with rejection-controlled EM
        rejection = (pcut>0)*ifelse((0.9^it.em)>=pcut,(0.9^it.em),pcut)
        DT[,c('Rejection1','Rejection2') := list(PostProb1,PostProb2)]
        DT[PostProb1<rejection,Rejection1 := rbinom(.N,1,prob=PostProb1/rejection)*rejection]
        DT[PostProb2<rejection,Rejection2 := rbinom(.N,1,prob=PostProb2/rejection)*rejection]

        ### Updating the vectorized dataset
        DTvec[,c('Rejection1','Rejection2') := .(rep(DT[,Rejection1],N),rep(DT[,Rejection2],N))]

        ### Aggregating data
        dt1 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(Rejection1>0)',agg = 'Rejection1')
        dt2 <- agg(data = DTvec,data.unique = DTvec.unique,rows = '(Rejection2>0)',agg = 'Rejection2')

        ### Calculating MLEs
        tryCatch({assign('model1',optim(par=inv.par(psi1.k,model='zinb'),fn=glm.zinb,gr = deriv.glm.zinb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi,rep(-Inf,ncolControl)),
                                        Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=F]),offset.vec=dt1[,offset],weights.vec=dt1[,weights]))},
                 error=function(e){model1<<-list();model1[['par']]<<-inv.par(psi1.k,model='zinb');model1[['convergence']]<<-99})
        tryCatch({assign('model2',optim(par=inv.par(psi2.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                        Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=F]),offset.vec=dt2[,offset],weights.vec=dt2[,weights]))},
                     error=function(e){model2<<-list();model2[['par']]<<-inv.par(psi2.k,model='nb');model2[['convergence']]<<-99})

        rm(dt1);rm(dt2)

        ### Saving parameters
        psi1.k1 = inv.par(model1$par,model='zinb');names(psi1.k1) = names(psi1.k)
        psi2.k1 = inv.par(model2$par,model='nb')
        psi.k1 = c(psi1.k1,psi2.k1)

        # Updating parameter history
        theta.k1 = c(pi.k1,gamma.k1,psi.k1)
        names(theta.k1) = names(theta.k)
        theta.k = theta.k1
        parlist[[it.em]] = c(it=it.em,Q=Q(as.matrix(DT[,.(PostProb1,PostProb2)]),as.matrix(DT[,.(JoinProb11,JoinProb12,JoinProb21,JoinProb22)]),loglik,pi.k1,gamma.k1),
                             error=error.em[1],theta.k1,m1=model1$convergence,m2=model2$convergence)

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

    if(!quiet){cat('\nDone!\n')}
    return(list('Pi'=pi.k1,'Gamma'=gamma.k1,'Psi'=psi.k1,'Prob'=DT[,.(PostProb1,PostProb2)],'LogF'=logF,'LogB'=logB,'Loglik'=loglik,'Parhist'=as.data.frame(do.call(rbind,parlist)),
                'Mean'=mu,'Viterbi'=zlist[[it.em]]))
}

ZIMHMM = function(ChIP,Control,offset,random,control)
{
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
        model.list = HMM.init(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
                              Control=NULL,
                              offset.init=DT[,rowMeans(.SD),.SDcols = paste0('offset.',1:N)],pcut=pcut)
    } else{
        model.list = HMM.init(ChIP.init=DT[,rowSums(.SD),.SDcols = paste0('ChIP.',1:N)],
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

ZIHMM.parallel = function(chromosome,ChIP,Control,offset,control){
    ChIP.sub = subset(ChIP,chr==chromosome);ChIP.sub$chr = NULL; ChIP.sub = as.matrix(ChIP.sub)
    Control.sub = subset(Control,chr==chromosome);Control.sub$chr = NULL; Control.sub = as.matrix(Control.sub)
    offset.sub = subset(offset,chr==chromosome);offset.sub$chr = NULL; offset.sub = as.matrix(offset.sub)
    tmp = ZIHMM(ChIP = ChIP.sub,Control = Control.sub,offset = offset.sub,control=control)
    return(tmp)
}

ZIMHMM.parallel = function(chromosome,ChIP,Control,offset,random,control){
    ChIP.sub = subset(ChIP,chr==chromosome);ChIP.sub$chr = NULL; ChIP.sub = as.matrix(ChIP.sub)
    Control.sub = subset(Control,chr==chromosome);Control.sub$chr = NULL; Control.sub = as.matrix(Control.sub)
    offset.sub = subset(offset,chr==chromosome);offset.sub$chr = NULL; offset.sub = as.matrix(offset.sub)
    tmp = ZIMHMM(ChIP = ChIP.sub,Control = Control.sub,offset = offset.sub,random = random,control = control)
    return(tmp)
}
