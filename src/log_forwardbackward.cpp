// Forward-backward algorithm for non-mixture hidden Markov models using log-space

#include "seqHMM.h"
// [[Rcpp::export]]

List log_forwardbackward(const arma::mat& transition_, const arma::cube& emission_, 
  const arma::vec& init_, const arma::ucube& obs, bool forwardonly, unsigned int threads) {

  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k

  log_internalForward(transition, emission, init, obs, alpha, threads);

  if (forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return List::create(Named("forward_probs") = wrap(alpha), Named("backward_probs") = wrap(beta));
  }

}
