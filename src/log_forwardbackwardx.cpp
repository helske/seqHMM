// Forward-backward algorithm for mixture hidden Markov models using log-space

#include "seqHMM.h"
// [[Rcpp::export]]

List log_forwardbackwardx(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, const arma::mat& coef, const arma::mat& X,
    const arma::ivec& numberOfStates, bool forwardonly, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);

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
