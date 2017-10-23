// log-likelihood of MHMM using log-space
#include <RcppArmadillo.h>
#include "reparma.h"

// [[Rcpp::export]]

Rcpp::NumericVector logLikMixHMM(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
  const arma::uvec& numberOfStates, unsigned int threads) {
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::wrap(-arma::datum::inf);
  }
  weights.each_row() /= sum(weights, 0);
  
  arma::vec ll(obs.n_slices);
  arma::sp_mat transition_t(transition.t());
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition_t, numberOfStates)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      arma::vec alpha = init % reparma(weights.col(k), numberOfStates);
      
      for (unsigned int r = 0; r < obs.n_rows; r++) {
        alpha %= emission.slice(r).col(obs(r, 0, k));
      }
      
      double tmp = sum(alpha);
      ll(k) = log(tmp);
      alpha /= tmp;
      
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        alpha = transition_t * alpha;
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          alpha %= emission.slice(r).col(obs(r, t, k));
        }
        
        tmp = sum(alpha);
        ll(k) += log(tmp);
        alpha /= tmp;
      }
    }
    return Rcpp::wrap(ll);
}
