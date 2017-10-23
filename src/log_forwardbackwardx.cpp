// Forward-backward algorithm for mixture hidden Markov models using log-space

#include "log_forward_backward.h"
#include "reparma.h"
// [[Rcpp::export]]

Rcpp::List log_forwardbackwardx(const arma::mat& transition_, 
  const arma::cube& emission_, const arma::vec& init_,
  const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
    const arma::uvec& numberOfStates, bool forwardonly, unsigned int threads) {

  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);
  
  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= arma::sum(weights, 0);
  weights = log(weights);

  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init + reparma(weights.col(k), numberOfStates);
  }
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  log_internalForwardx(transition, emission, initk, obs, alpha, threads);

  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta));
  }

  return Rcpp::wrap(alpha);
}
