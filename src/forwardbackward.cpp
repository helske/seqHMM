#include "seqHMM.h"

// [[Rcpp::export]]

List forwardbackward(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, bool forwardonly, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_rows); //n,k

  internalForward(transition, emission, init, obs, alpha, scales, threads);

  if (forwardonly) {
    return List::create(Named("forward_probs") = wrap(alpha),
        Named("scaling_factors") = wrap(scales));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k
    internalBackward(transition, emission, obs, beta, scales, threads);
    return List::create(Named("forward_probs") = wrap(alpha), Named("backward_probs") = wrap(beta),
        Named("scaling_factors") = wrap(scales));
  }

}
