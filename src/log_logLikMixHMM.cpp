#include "seqHMM.h"
// [[Rcpp::export]]

NumericVector log_logLikMixHMM(NumericVector transitionMatrix, NumericVector emissionArray,
  NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, NumericMatrix X_,
  IntegerVector numberOfStates, int threads) {
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(), q, coefs.ncol(), false);
  coef.col(0).zeros();
  arma::mat X(X_.begin(), obs.n_rows, q, false);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    warning(
      "Coefficients of covariates resulted non-finite cluster probabilities. Returning -Inf.");
    return wrap(-arma::math::inf());
    
  }
  weights.each_row() /= sum(weights, 0);
  
  weights = log(weights);
  transition = log(transition);
  emission = log(emission);
  init = log(init);
  
  NumericVector ll(obs.n_rows);
  
#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition, numberOfStates)
    for (int k = 0; k < obs.n_rows; k++) {
      arma::vec alpha = init + reparma(weights.col(k), numberOfStates);
      for (int r = 0; r < obs.n_slices; r++) {
        alpha += emission.slice(r).col(obs(k, 0, r));
      }
      arma::vec alphatmp(emission.n_rows);
      for (int t = 1; t < obs.n_cols; t++) {
        for (int i = 0; i < emission.n_rows; i++) {
          alphatmp(i) = logSumExp(alpha + transition.col(i));
          for (int r = 0; r < obs.n_slices; r++) {
            alphatmp(i) += emission(i, obs(k, t, r), r);
          }
        }
        alpha = alphatmp;
      }
      ll(k) = logSumExp(alpha);
    }
    
    return ll;
}

