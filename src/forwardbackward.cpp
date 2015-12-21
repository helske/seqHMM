// Forward-backward algorithm for non-mixture hidden Markov models
#include "seqHMM.h"

// [[Rcpp::export]]

List forwardbackward(const arma::mat& transition, NumericVector emissionArray,
    const arma::vec& init, IntegerVector obsArray, bool forwardonly, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false, true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //n,k

  internalForward(transition, emission, init, obs, alpha, scales, threads);
  

  if(!scales.is_finite()) {
    Rcpp::stop("Scaling factors contain non-finite values. \n Check the model or try using the log-space version of the algorithm.");
  }
  double min_sf = scales.min();
  if (min_sf < 1e-150) {
    Rcpp::warning("Smallest scaling factor was %e, results can be numerically unstable.", min_sf);
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
