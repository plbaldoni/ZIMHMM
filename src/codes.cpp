// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <RcppEigen.h>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

//[[Rcpp::export]]
NumericVector Viterbi(NumericMatrix LOGF, NumericVector P, NumericMatrix GAMMA){
    int M = LOGF.nrow();
    int K = LOGF.ncol();

    NumericMatrix LOGV(M,K);
    NumericMatrix V(M,K);

    LOGV(0,_) = log(P) + LOGF(0,_);
    V(0,_) = P*exp(LOGF(0,_));
    V(0,_) = V(0,_)/sum(V(0,_));

    NumericMatrix COND(K,K);
    NumericMatrix MAXCOND(K);

    for(int j = 0; j < (M-1); j++){
        for(int k = 0; k < K; k++){
            COND(k,_) = V(j,_)*GAMMA(_,k);
            MAXCOND[k] = max(COND(k,_));
        }
        V(j+1,_) = MAXCOND*exp(LOGF(j+1,_));
        V(j+1,_) = V(j+1,_)/sum(V(j+1,_));
    }

    NumericVector S(M);
    NumericVector AUXS(K);
    S[M-1] = which_max(V(M-1,_));

    for(int j = (M-2); j >= 0; j--){
        AUXS = V(j,_)*GAMMA(_,S[j+1]);
        S[j] = which_max(AUXS);
    }

    return(S);
}

//[[Rcpp::export]]
NumericVector hmm_logF(NumericVector logf1, NumericVector logf2, NumericVector pi, NumericMatrix gamma){
    //f1 = background log-density evaluated at response vector
    //f2 = enrichment log-density evaluated at response vector
    //pi = probability vector
    //gamma = transition probability

    int N = logf1.size(); //Number of observations
    int K = pi.size(); //Number of states

    NumericMatrix logF(N,K); //Matrix f(i,k) (Equation A.1. JASA paper)
    NumericMatrix F(N,K); //Matrix f(i,k) (Equation A.1. JASA paper)

    //Initializing logF
    logF(0,0) = log(pi[0]) + logf1[0];
    logF(0,1) = log(pi[1]) + logf2[0];

    //Forward loop
    for(int i = 1; i < N; i++){
        NumericVector a(2);
        a(0) = logF(i-1,0) + log(gamma(0,0)) + logf1[i];
        a(1) = logF(i-1,1) + log(gamma(1,0)) + logf1[i];
        double maxa = max(a);

        NumericVector c(2);
        c(0) = logF(i-1,0) + log(gamma(0,1)) + logf2[i];
        c(1) = logF(i-1,1) + log(gamma(1,1)) + logf2[i];
        double maxc = max(c);

        logF(i,0) = maxa + log(exp(a(0)-maxa)+exp(a(1)-maxa));
        logF(i,1) = maxc + log(exp(c(0)-maxc)+exp(c(1)-maxc));
    }
    return logF;
}

//[[Rcpp::export]]
NumericVector hmm_logB(NumericVector logf1, NumericVector logf2, NumericVector pi, NumericMatrix gamma){
    //f1 = background log-density evaluated at response vector
    //f2 = enrichment log-density evaluated at response vector
    //pi = probability vector
    //gamma = transition probability

    int N = logf1.size(); //Number of observations
    int K = pi.size(); //Number of states

    NumericMatrix logB(N,K); //Matrix b(i,k) (Equation A.2. JASA paper)
    NumericMatrix B(N,K); //Matrix b(i,k) (Equation A.2. JASA paper)

    //Initializing logB
    logB(N-1,0) = log(1);
    logB(N-1,1) = log(1);

    //Backward loop
    for(int j = (N-2); j >= 0; j--){
        NumericVector a(2);
        a(0) = logB(j+1,0) + log(gamma(0,0)) + logf1[j+1];
        a(1) = logB(j+1,1) + log(gamma(0,1)) + logf2[j+1];
        double maxa = max(a);

        NumericVector c(2);
        c(0) = logB(j+1,0) + log(gamma(1,0)) + logf1[j+1];
        c(1) = logB(j+1,1) + log(gamma(1,1)) + logf2[j+1];
        double maxc = max(c);

        logB(j,0) = maxa + log(exp(a(0)-maxa) + exp(a(1)-maxa));
        logB(j,1) = maxc + log(exp(c(0)-maxc) + exp(c(1)-maxc));
    }
    return logB;
}

//[[Rcpp::export]]
NumericVector hmm_P1(NumericMatrix logF, NumericMatrix logB){
    int N = logF.nrow(); //Number of observations
    int K = logF.ncol(); //Number of states

    NumericMatrix logP(N,K); //Final probability matrix (output)
    NumericMatrix P(N,K); //Final probability matrix (output)

    NumericVector d(2);
    d(0) = logF(N-1,0);
    d(1) = logF(N-1,1);
    double maxd = max(d);

    //Output loop
    for(int i = 0; i < N;i++){
        for(int j = 0; j < K; j++){
            logP(i,j) = logF(i,j) + logB(i,j) - maxd - log(exp(d(0)-maxd) + exp(d(1)-maxd));
            P(i,j) = exp(logP(i,j));
        }
    }
    return P;
}

//[[Rcpp::export]]
NumericVector hmm_P2(NumericMatrix logF, NumericMatrix logB,NumericVector logf1, NumericVector logf2, NumericMatrix gamma){
    //f1 = background log-density evaluated at response vector
    //f2 = enrichment log-density evaluated at response vector
    //pi = probability vector
    //gamma = transition probability

    int N = logF.nrow(); //Number of observations
    int K = logF.ncol(); //Number of states

    NumericMatrix logP(N,K*K); //Final probability matrix (output)
    NumericMatrix P(N,K*K); //Final probability matrix (output)

    NumericVector d(2);
    d(0) = logF(N-1,0);
    d(1) = logF(N-1,1);
    double maxd = max(d);

    //Output loop
    for(int i = 0; i < N;i++){
        if(i >= 1){
            logP(i,0) = logF(i-1,0) + log(gamma(0,0)) + logf1[i] + logB(i,0) - maxd - log(exp(d(0)-maxd) + exp(d(1)-maxd));
            logP(i,1) = logF(i-1,0) + log(gamma(0,1)) + logf2[i] + logB(i,1) - maxd - log(exp(d(0)-maxd) + exp(d(1)-maxd));
            logP(i,2) = logF(i-1,1) + log(gamma(1,0)) + logf1[i] + logB(i,0) - maxd - log(exp(d(0)-maxd) + exp(d(1)-maxd));
            logP(i,3) = logF(i-1,1) + log(gamma(1,1)) + logf2[i] + logB(i,1) - maxd - log(exp(d(0)-maxd) + exp(d(1)-maxd));

            P(i,0) = exp(logP(i,0));
            P(i,1) = exp(logP(i,1));
            P(i,2) = exp(logP(i,2));
            P(i,3) = exp(logP(i,3));
        }
    }
    //1st col = P(1,1), 2nd col = P(1,2), 3rd col = P(2,1), 4th col = P(2,2)
    return P;
}

//[[Rcpp::export]]
arma::mat my_dnbinom(arma::mat x, arma::mat mu, double size, int lg){
    int M = x.size();
    arma::mat out(M,1);
    for(int i = 0; i < M; i++){
        out(i,0) = R::dnbinom_mu(x[i],size,mu[i],lg);
    }
    return(out);
}

//[[Rcpp::export]]
arma::mat MHMMmean(arma::mat XMAT, arma::mat BETA, arma::mat RANDOM, arma::mat OFFSETVEC, int N, int M, int K)
{
    arma::mat MU(M,K*N); //Final mean matrix

    int index = 0;
    for(int i = 0; i < N; i++){
        int from = i*M;
        int to = from + (M - 1);
        arma::mat subXMAT = XMAT.rows(from,to);
        arma::mat subRANDOM = RANDOM.rows(from,to);
        arma::mat subOFFSETVEC = OFFSETVEC.rows(from,to);
        for(int k = 0; k < K; k++){
            arma::mat subBETA = BETA.row(k);
            arma::mat subMU = exp(subXMAT*subBETA.t()+subRANDOM+subOFFSETVEC);
            MU.col(index) = subMU;
            index++;
        }
    }
    return MU;
}

//[[Rcpp::export]]
arma::mat MHMMLik(arma::vec YVEC, arma::mat ZEROINFL, arma::mat MU, arma::vec DISP,int N, int M, int K){
    arma::mat LL(M,K); LL.zeros();

    //NumericVector onevec(M,1.0);
    //NumericVector zerovec(M,0.0);
    //NumericMatrix wrap_YVEC = wrap(YVEC);

    //arma::mat YVEC0 = ifelse(wrap_YVEC==zerovec,onevec,zerovec);
    arma::mat YVEC0(N*M,1); YVEC0.zeros();
    arma::uvec ids = find(YVEC == 0);
    YVEC0.rows(ids).fill(1);
    arma::mat ZEROINFLvec(M*N,1); ZEROINFLvec.col(0) = vectorise(ZEROINFL);

    //Background
    NumericVector idxNV = wrap(seq(0,(N-1))); idxNV = idxNV*K;
    arma::uvec idx = as<arma::uvec>(idxNV);
    arma::mat muidx(M*N,1); muidx.col(0) = vectorise(MU.cols(idx));
    arma::mat zeroidx(N*M,1); zeroidx.zeros();
    arma::mat LogLik00(N*M,1); LogLik00.col(0) = my_dnbinom(zeroidx,muidx,DISP[0],1);
    arma::mat LogLik01(N*M,1); LogLik01.col(0) = my_dnbinom(YVEC,muidx,DISP[0],1);

    arma::mat LLvecaux(N*M,1); LLvecaux.col(0) = YVEC0%log(ZEROINFLvec + exp(log(1-ZEROINFLvec)+LogLik00)) + (1-YVEC0)%(log(1-ZEROINFLvec) + LogLik01);
    arma::mat LLmataux = reshape(LLvecaux,M,N);
    NumericMatrix LLmataux1 = wrap(LLmataux);
    arma::mat LLmataux2 = rowSums(LLmataux1);
    LL.col(0) = LLmataux2;

    //Enrichment
    idxNV = idxNV + 1;
    idx = as<arma::uvec>(idxNV);
    muidx.col(0) = vectorise(MU.cols(idx));
    arma::mat LogLik1(N*M,1); LogLik1.col(0) = my_dnbinom(YVEC,muidx,DISP[1],1);
    LLmataux = reshape(LogLik1,M,N);
    LLmataux1 = wrap(LLmataux);
    LLmataux2 = rowSums(LLmataux1);
    LL.col(1) = LLmataux2;

    return(LL);
}

// [[Rcpp::export]]
double integrand(NumericVector U, arma::vec YVEC, arma::mat XMAT, arma::mat BETA, arma::vec DISP,
                 NumericVector P, NumericMatrix GAMMA, arma::vec OFFSETVEC,
                 arma::mat ZEROINFL,NumericVector W, double SIGMA2){

    int N = U.size();
    int M = YVEC.size()/N;
    int K = P.size();

    NumericVector RANDOMrcpp = sqrt(SIGMA2)*rep_each(U,M)*W;
    arma::mat RANDOM(N*M,1); RANDOM = as<arma::vec>(RANDOMrcpp);

    arma::mat MU = MHMMmean(XMAT,BETA,RANDOM,OFFSETVEC,N,M,K);
    arma::mat LogLik = MHMMLik(YVEC,ZEROINFL,MU,DISP,N,M,K);

    NumericVector logFj1(K);
    NumericVector logFj2(K);
    NumericMatrix Bj1(K,K);
    NumericMatrix Bj2(K,K);
    NumericMatrix Bj(K,K);

    for(int j = 2; j < M; j++){
        if(j==2){
            logFj1(0) = LogLik(j-1,0);//R::dnbinom_mu(Y[j-1],PHI[0],exp(B0[0]+B1[0]*X[j-1]+U),1);
            logFj1(1) = LogLik(j-1,1);//R::dnbinom_mu(Y[j-1],PHI[1],exp(B0[1]+B1[1]*X[j-1]+U),1);
            Bj1(0,0) = log(GAMMA(0,0)) + logFj1[0];Bj1(0,1) = log(GAMMA(0,1)) + logFj1[1];
            Bj1(1,0) = log(GAMMA(1,0)) + logFj1[0];Bj1(1,1) = log(GAMMA(1,1)) + logFj1[1];
        } else{
            arma::mat Bj1 = log(as<arma::mat>(clone(Bj)));
        }

        logFj2(0) = LogLik(j,0);//R::dnbinom_mu(Y[j],PHI[0],exp(B0[0]+B1[0]*X[j]+U),1);
        logFj2(1) = LogLik(j,1);//R::dnbinom_mu(Y[j],PHI[1],exp(B0[1]+B1[1]*X[j]+U),1);
        Bj2(0,0) = log(GAMMA(0,0)) + logFj2[0];Bj2(0,1) = log(GAMMA(0,1)) + logFj2[1];
        Bj2(1,0) = log(GAMMA(1,0)) + logFj2[0];Bj2(1,1) = log(GAMMA(1,1)) + logFj2[1];

        Bj(0,0) = exp(Bj1(0,0)+Bj2(0,0)) + exp(Bj1(0,1)+Bj2(1,0));
        Bj(0,1) = exp(Bj1(0,0)+Bj2(0,1)) + exp(Bj1(0,1)+Bj2(1,1));
        Bj(1,0) = exp(Bj1(1,0)+Bj2(0,0)) + exp(Bj1(1,1)+Bj2(1,0));
        Bj(1,1) = exp(Bj1(1,0)+Bj2(0,1)) + exp(Bj1(1,1)+Bj2(1,1));
    }
    NumericVector A(2);
    A[0] = log(P[0]) + LogLik(0,0);//R::dnbinom_mu(Y[0],PHI[0],exp(B0[0]+B1[0]*X[0]+U),1);
    A[1] = log(P[1]) + LogLik(0,1);//R::dnbinom_mu(Y[0],PHI[1],exp(B0[1]+B1[1]*X[0]+U),1);

    NumericVector AB(2);
    AB[0] = exp(A[0]+log(Bj(0,0)))+exp(A[1]+log(Bj(1,0)));
    AB[1] = exp(A[0]+log(Bj(0,1)))+exp(A[1]+log(Bj(1,1)));

    double result = log(sum(AB)) + sum(dnorm(U,0.0,1.0,1));
    return(result);
}

//[[Rcpp::export]]
NumericVector generateHMM(NumericMatrix GAMMA,int m){
    //generateHMM simulates a Markov Chain given a matrix of transition probabilities
    //This function simulates a sequence of hidden states Z based on the transition probabilities GAMMA
    //This function is used when simulating data
    NumericVector Z(m,0.0);
    Z[0] = 1;
    for(int i = 1; i < m; i ++)
    {
        double rand = R::runif(0,1);
        double sump = 0;
        int flag = 0;
        int k = 0;
        while(flag == 0){
            sump = sump + GAMMA((Z[i-1]-1),k);
            if(rand<=sump){
                Z[i] = (k+1);
                flag++;
            }
            k++;
        }
    }
    return Z;
}

