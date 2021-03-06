// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Viterbi
NumericVector Viterbi(NumericMatrix LOGF, NumericVector P, NumericMatrix GAMMA);
RcppExport SEXP _ZIMHMM_Viterbi(SEXP LOGFSEXP, SEXP PSEXP, SEXP GAMMASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type LOGF(LOGFSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type GAMMA(GAMMASEXP);
    rcpp_result_gen = Rcpp::wrap(Viterbi(LOGF, P, GAMMA));
    return rcpp_result_gen;
END_RCPP
}
// hmm_logF
NumericVector hmm_logF(NumericVector logf1, NumericVector logf2, NumericVector pi, NumericMatrix gamma);
RcppExport SEXP _ZIMHMM_hmm_logF(SEXP logf1SEXP, SEXP logf2SEXP, SEXP piSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logf1(logf1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type logf2(logf2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_logF(logf1, logf2, pi, gamma));
    return rcpp_result_gen;
END_RCPP
}
// hmm_logB
NumericVector hmm_logB(NumericVector logf1, NumericVector logf2, NumericVector pi, NumericMatrix gamma);
RcppExport SEXP _ZIMHMM_hmm_logB(SEXP logf1SEXP, SEXP logf2SEXP, SEXP piSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logf1(logf1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type logf2(logf2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_logB(logf1, logf2, pi, gamma));
    return rcpp_result_gen;
END_RCPP
}
// hmm_P1
NumericVector hmm_P1(NumericMatrix logF, NumericMatrix logB);
RcppExport SEXP _ZIMHMM_hmm_P1(SEXP logFSEXP, SEXP logBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type logF(logFSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type logB(logBSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_P1(logF, logB));
    return rcpp_result_gen;
END_RCPP
}
// hmm_P2
NumericVector hmm_P2(NumericMatrix logF, NumericMatrix logB, NumericVector logf1, NumericVector logf2, NumericMatrix gamma);
RcppExport SEXP _ZIMHMM_hmm_P2(SEXP logFSEXP, SEXP logBSEXP, SEXP logf1SEXP, SEXP logf2SEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type logF(logFSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type logB(logBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type logf1(logf1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type logf2(logf2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_P2(logF, logB, logf1, logf2, gamma));
    return rcpp_result_gen;
END_RCPP
}
// my_dnbinom
arma::mat my_dnbinom(arma::mat x, arma::mat mu, double size, int lg);
RcppExport SEXP _ZIMHMM_my_dnbinom(SEXP xSEXP, SEXP muSEXP, SEXP sizeSEXP, SEXP lgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type lg(lgSEXP);
    rcpp_result_gen = Rcpp::wrap(my_dnbinom(x, mu, size, lg));
    return rcpp_result_gen;
END_RCPP
}
// MHMMmean
arma::mat MHMMmean(arma::mat XMAT, arma::mat BETA, arma::mat RANDOM, arma::mat OFFSETVEC, int N, int M, int K);
RcppExport SEXP _ZIMHMM_MHMMmean(SEXP XMATSEXP, SEXP BETASEXP, SEXP RANDOMSEXP, SEXP OFFSETVECSEXP, SEXP NSEXP, SEXP MSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XMAT(XMATSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RANDOM(RANDOMSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type OFFSETVEC(OFFSETVECSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(MHMMmean(XMAT, BETA, RANDOM, OFFSETVEC, N, M, K));
    return rcpp_result_gen;
END_RCPP
}
// MHMMLik
arma::mat MHMMLik(arma::vec YVEC, arma::mat ZEROINFL, arma::mat MU, arma::vec DISP, int N, int M, int K);
RcppExport SEXP _ZIMHMM_MHMMLik(SEXP YVECSEXP, SEXP ZEROINFLSEXP, SEXP MUSEXP, SEXP DISPSEXP, SEXP NSEXP, SEXP MSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type YVEC(YVECSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ZEROINFL(ZEROINFLSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type MU(MUSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type DISP(DISPSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(MHMMLik(YVEC, ZEROINFL, MU, DISP, N, M, K));
    return rcpp_result_gen;
END_RCPP
}
// integrand
double integrand(NumericVector U, arma::vec YVEC, arma::mat XMAT, arma::mat BETA, arma::vec DISP, NumericVector P, NumericMatrix GAMMA, arma::vec OFFSETVEC, arma::mat ZEROINFL, NumericVector W, double SIGMA2);
RcppExport SEXP _ZIMHMM_integrand(SEXP USEXP, SEXP YVECSEXP, SEXP XMATSEXP, SEXP BETASEXP, SEXP DISPSEXP, SEXP PSEXP, SEXP GAMMASEXP, SEXP OFFSETVECSEXP, SEXP ZEROINFLSEXP, SEXP WSEXP, SEXP SIGMA2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YVEC(YVECSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XMAT(XMATSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type DISP(DISPSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type GAMMA(GAMMASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type OFFSETVEC(OFFSETVECSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ZEROINFL(ZEROINFLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type SIGMA2(SIGMA2SEXP);
    rcpp_result_gen = Rcpp::wrap(integrand(U, YVEC, XMAT, BETA, DISP, P, GAMMA, OFFSETVEC, ZEROINFL, W, SIGMA2));
    return rcpp_result_gen;
END_RCPP
}
// generateHMM
NumericVector generateHMM(NumericMatrix GAMMA, int m);
RcppExport SEXP _ZIMHMM_generateHMM(SEXP GAMMASEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type GAMMA(GAMMASEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(generateHMM(GAMMA, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ZIMHMM_Viterbi", (DL_FUNC) &_ZIMHMM_Viterbi, 3},
    {"_ZIMHMM_hmm_logF", (DL_FUNC) &_ZIMHMM_hmm_logF, 4},
    {"_ZIMHMM_hmm_logB", (DL_FUNC) &_ZIMHMM_hmm_logB, 4},
    {"_ZIMHMM_hmm_P1", (DL_FUNC) &_ZIMHMM_hmm_P1, 2},
    {"_ZIMHMM_hmm_P2", (DL_FUNC) &_ZIMHMM_hmm_P2, 5},
    {"_ZIMHMM_my_dnbinom", (DL_FUNC) &_ZIMHMM_my_dnbinom, 4},
    {"_ZIMHMM_MHMMmean", (DL_FUNC) &_ZIMHMM_MHMMmean, 7},
    {"_ZIMHMM_MHMMLik", (DL_FUNC) &_ZIMHMM_MHMMLik, 7},
    {"_ZIMHMM_integrand", (DL_FUNC) &_ZIMHMM_integrand, 11},
    {"_ZIMHMM_generateHMM", (DL_FUNC) &_ZIMHMM_generateHMM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ZIMHMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
