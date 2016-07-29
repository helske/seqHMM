// Forward-backward algorithm for mixture hidden Markov models using log-space

#include "seqHMM.h"
// [[Rcpp::export]]

List log_forwardbackwardx(arma::mat transition, arma::cube emission, arma::vec init,
  const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
    const arma::uvec& numberOfStates, bool forwardonly, unsigned int threads) {

  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= arma::sum(weights, 0);
  weights = log(weights);
  transition = log(transition);
  emission = log(emission);
  init = log(init);

  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init + reparma(weights.col(k), numberOfStates);
  }
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  log_internalForwardx(transition, emission, initk, obs, alpha, threads);

  if (forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return List::create(Named("forward_probs") = wrap(alpha), Named("backward_probs") = wrap(beta));
  }

  return wrap(alpha);
}
