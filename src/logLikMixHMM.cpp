// log-likelihood of MHMM using log-space
#include "seqHMM.h"

// [[Rcpp::export]]

NumericVector logLikMixHMM(const arma::mat& transition, NumericVector emissionArray,
  const arma::vec& init, IntegerVector obsArray, const arma::mat& coef, const arma::mat& X,
  const arma::ivec& numberOfStates, int threads) {
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false, true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return wrap(-arma::math::inf());
  }
  weights.each_row() /= sum(weights, 0);
  
  arma::vec ll(obs.n_slices);
  arma::sp_mat transition_t(transition.t());
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition_t, numberOfStates)
    for (int k = 0; k < obs.n_slices; k++) {
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
    return wrap(ll);
}
