// log-likelihood of MHMM using log-space

#include "seqHMM.h"
// [[Rcpp::export]]

NumericVector log_logLikMixHMM(NumericVector transitionMatrix, NumericVector emissionArray,
  NumericVector initialProbs, IntegerVector obsArray, const arma::mat& coef, const arma::mat& X,
  const arma::ivec& numberOfStates, int threads) {
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return wrap(-arma::math::inf());
  }
  weights.each_row() /= sum(weights, 0);
  
  weights = log(weights);
  transition = log(transition);
  emission = log(emission);
  init = log(init);
  
  arma::vec ll(obs.n_slices);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition, numberOfStates)
    for (int k = 0; k < obs.n_slices; k++) {
      arma::vec alpha = init + reparma(weights.col(k), numberOfStates);
      for (int r = 0; r < obs.n_rows; r++) {
        alpha += emission.slice(r).col(obs(r, 0, k));
      }
      arma::vec alphatmp(emission.n_rows);
      for (int t = 1; t < obs.n_cols; t++) {
        for (int i = 0; i < emission.n_rows; i++) {
          alphatmp(i) = logSumExp(alpha + transition.col(i));
          for (int r = 0; r < obs.n_rows; r++) {
            alphatmp(i) += emission(i, obs(r, t, k), r);
          }
        }
        alpha = alphatmp;
      }
      ll(k) = logSumExp(alpha);
    }
    
    return wrap(ll);
}

