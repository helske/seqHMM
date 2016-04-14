//stand-alone coefficient estimation
#include "seqHMM.h"

// [[Rcpp::export]]
List estimate_coefs(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
         IntegerVector obsArray, const arma::ivec& nSymbols, NumericMatrix coefs, const arma::mat& X,
         const arma::ivec& numberOfStates, int itermax, double tol, int trace, int threads) {
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);
  
  arma::mat coef(coefs.begin(), coefs.nrow(), coefs.ncol());
  coef.col(0).zeros();
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return List::create(Named("error") = 3);
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
    return List::create(Named("error") = 1);
  }
  internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
  if(!beta.is_finite()) {
    return List::create(Named("error") = 2);
  }
  double min_sf = scales.min();
  if (min_sf < 1e-150) {
    Rcpp::warning("Smallest scaling factor was %e, results can be numerically unstable. ", min_sf);
  }
  
  double sumlogLik = arma::accu(log(scales));
  
  if (trace > 0) {
    Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
  }
  
  double change = tol + 1.0;
  int iter = 0;
  
  arma::ivec cumsumstate = arma::cumsum(numberOfStates);
  
  while ((change > tol) & (iter < itermax)) {
    iter++;
  
  unsigned int error = optCoef(weights, obs, emission, initk, beta, scales, coef, X, cumsumstate,
                                   numberOfStates, trace);
  if (error != 0) {
    return List::create(Named("error") = error);
  }
  
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }
  
  internalForwardx(sp_trans.t(), emission, initk, obs, alpha, scales, threads);
  if(!scales.is_finite()) {
    return List::create(Named("error") = 1);
  }
  internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
  if(!beta.is_finite()) {
    return List::create(Named("error") = 2);
  }
  double min_sf = scales.min();
  if (min_sf < 1e-150) {
    Rcpp::warning("Smallest scaling factor was %e, results can be numerically unstable. ", min_sf);
  }
  double tmp = arma::accu(log(scales));
  change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
  sumlogLik = tmp;
  if (trace > 1) {
    Rcout << "iter: " << iter;
    Rcout << " logLik: " << sumlogLik;
    Rcout << " relative change: " << change << std::endl;
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
  return List::create(Named("coefficients") = wrap(coef), Named("error") = 0);
}
