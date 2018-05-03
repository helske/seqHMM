// log-likelihood of MHMM using log-space

#include "logsumexp.h"
#include "reparma.h"
// [[Rcpp::export]]

Rcpp::NumericVector log_logLikMixHMM(arma::mat transition, arma::cube emission, arma::vec init,
  const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
  const arma::uvec& numberOfStates, unsigned int threads) {
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::wrap(-arma::datum::inf);
  }
  weights.each_row() /= sum(weights, 0);
  
  weights = log(weights);
  transition = log(transition);
  emission = log(emission);
  init = log(init);
  
  arma::vec ll(obs.n_slices);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition, numberOfStates)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      arma::vec alpha = init + reparma(weights.col(k), numberOfStates);
      for (unsigned int r = 0; r < obs.n_rows; r++) {
        alpha += emission.slice(r).col(obs(r, 0, k));
      }
      arma::vec alphatmp(emission.n_rows);
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        for (unsigned int i = 0; i < emission.n_rows; i++) {
          alphatmp(i) = logSumExp(alpha + transition.col(i));
          for (unsigned int r = 0; r < obs.n_rows; r++) {
            alphatmp(i) += emission(i, obs(r, t, k), r);
          }
        }
        alpha = alphatmp;
      }
      ll(k) = logSumExp(alpha);
    }
    
    return Rcpp::wrap(ll);
}

