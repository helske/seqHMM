// Forward-backward algorithm for non-mixture hidden Markov models
#include "seqHMM.h"

// [[Rcpp::export]]

List forwardbackward(const arma::mat& transition, const arma::cube& emission,
    const arma::vec& init, const arma::ucube& obs, bool forwardonly, unsigned int threads) {

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //n,k

  internalForward(transition, emission, init, obs, alpha, scales, threads);

  if(!scales.is_finite()) {
    Rcpp::stop("Scaling factors contain non-finite values. \n Check the model or try using the log-space version of the algorithm.");
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }
  
  if (forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha),
        Named("scaling_factors") = wrap(scales));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    internalBackward(transition, emission, obs, beta, scales, threads);
    if(!beta.is_finite()) {
      Rcpp::stop("Backward probabilities contain non-finite values. Check the model or try using the log-space version of the algorithm.");
    }
    return List::create(Named("forward_probs") = wrap(alpha), Named("backward_probs") = wrap(beta),
        Named("scaling_factors") = wrap(scales));
  }

}
