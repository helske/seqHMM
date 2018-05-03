// Forward-backward algorithm for non-mixture hidden Markov models using log-space

#include "log_forward_backward.h"
// [[Rcpp::export]]

Rcpp::List log_forwardbackward(const arma::mat& transition_, const arma::cube& emission_, 
  const arma::vec& init_, const arma::ucube& obs, bool forwardonly, unsigned int threads) {

  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k

  log_internalForward(transition, emission, init, obs, alpha, threads);

  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta));
  }

}
