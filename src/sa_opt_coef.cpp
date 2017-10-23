//stand-alone coefficient estimation
#include "optcoef.h"
#include "forward_backward.h"
#include "reparma.h"

// [[Rcpp::export]]
Rcpp::List estimate_coefs(const arma::mat& transition, const arma::cube& emission, 
  const arma::vec& init, const arma::ucube& obs, const arma::uvec& nSymbols, 
  arma::mat coef, const arma::mat& X, const arma::uvec& numberOfStates, 
  int itermax, double tol, int trace, unsigned int threads) {

  coef.col(0).zeros();
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 3);
  }
  weights.each_row() /= sum(weights, 0);

  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //m,n,k

  arma::sp_mat sp_trans(transition);
  internalForwardx(sp_trans.t(), emission, initk, obs, alpha, scales, threads);
  if(!scales.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 1);
  }
  internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
  if(!beta.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 2);
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }

  double sumlogLik = -arma::accu(log(scales));

  if (trace > 0) {
    Rcpp:: Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
  }

  double change = tol + 1.0;
  int iter = 0;

  arma::uvec cumsumstate = arma::cumsum(numberOfStates);

  while ((change > tol) & (iter < itermax)) {
    iter++;
  //
  // unsigned int error = optCoef(weights, obs, emission, initk, beta, scales, coef, X, cumsumstate,
  //                                  numberOfStates, trace);
  // if (error != 0) {
  //   return Rcpp::List::create(Rcpp::Named("error") = error);
  // }

  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }

  internalForwardx(sp_trans.t(), emission, initk, obs, alpha, scales, threads);
  if(!scales.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 1);
  }
  internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
  if(!beta.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 2);
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }
  double tmp = -arma::accu(log(scales));
  change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
  sumlogLik = tmp;
  if (trace > 1) {
    Rcpp::Rcout << "iter: " << iter;
    Rcpp::Rcout << " logLik: " << sumlogLik;
    Rcpp::Rcout << " relative change: " << change << std::endl;
  }

  }
  if (trace > 0) {
    if (iter == itermax) {
      Rcpp::Rcout << "EM algorithm stopped after reaching the maximum number of " << iter
                  << " iterations." << std::endl;
    } else {
      Rcpp::Rcout << "EM algorithm stopped after reaching the relative change of " << change;
      Rcpp::Rcout << " after " << iter << " iterations." << std::endl;
    }
    Rcpp::Rcout << "Final log-likelihood: " << sumlogLik << std::endl;

  }
  return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::wrap(coef), Rcpp::Named("error") = 0);
}
