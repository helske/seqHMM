// Forward-backward algorithm for non-mixture hidden Markov models using log-space

#include "seqHMM.h"
// [[Rcpp::export]]

List log_forwardbackward(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, bool forwardonly, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);

  init = log(init);
  transition = log(transition);
  emission = log(emission);

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
