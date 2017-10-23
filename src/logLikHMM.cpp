// log-likelihood of HMM
#include <RcppArmadillo.h>
// [[Rcpp::export]]
Rcpp::NumericVector logLikHMM(const arma::mat& transition, const arma::cube& emission,
    const arma::vec& init, const arma::ucube& obs, unsigned int threads) {

  arma::vec ll(obs.n_slices);
  arma::mat transition_t(transition.t());
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition_t)
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    arma::vec alpha = init;
    
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

