// Forward-backward algorithm for mixture hidden Markov models
#include "forward_backward.h"
#include "reparma.h"

// [[Rcpp::export]]

Rcpp::List forwardbackwardx(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube obs, const arma::mat& coef, const arma::mat& X,
  const arma::uvec& numberOfStates, bool forwardonly, unsigned int threads) {
  
  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= arma::sum(weights, 0);
  
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //m,n,k
  
  arma::sp_mat sp_trans(transition);
  internalForwardx(sp_trans.t(), emission, initk, obs, alpha, scales, threads);

  if(!scales.is_finite()) {
    Rcpp::stop("Scaling factors contain non-finite values. Check the model or try using the log-space version of the algorithm.");
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }
  
  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha),
      Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
    if(!beta.is_finite()) {
      Rcpp::stop("Backward probabilities contain non-finite values. Check the model or try using the log-space version of the algorithm.");
    }
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta),
      Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  }
  
  return Rcpp::wrap(alpha);
}
