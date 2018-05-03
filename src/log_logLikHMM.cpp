// log-likelihood of HMM using log-space

#include "logsumexp.h"
// [[Rcpp::export]]

Rcpp::NumericVector log_logLikHMM(const arma::mat& transition_, const arma::cube& emission_, 
  const arma::vec& init_, const arma::ucube& obs, unsigned int threads) {
  
  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);
  
  arma::vec ll(obs.n_slices);
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      arma::vec alpha = init;
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

