#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector logLikMixHMM(NumericVector transitionMatrix, NumericVector emissionArray,
  NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, NumericMatrix X_,
  IntegerVector numberOfStates, int threads = 1) {
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(), q, coefs.ncol(), false);
  coef.col(0).zeros();
  arma::mat X(X_.begin(), obs.n_rows, q, false);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    warning(
      "Coefficients of covariates resulted non-finite cluster probabilities. Returning -Inf.");
    return wrap(-arma::math::inf());
    
  }
  weights.each_row() /= sum(weights, 0);
  
  NumericVector ll(obs.n_rows);
  
#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition, numberOfStates)
    for (int k = 0; k < obs.n_rows; k++) {
      arma::vec alpha = init % reparma(weights.col(k), numberOfStates);
      
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        alpha %= emission.slice(r).col(obs(k, 0, r));
      }
      
      double tmp = sum(alpha);
      ll(k) = log(tmp);
      alpha /= tmp;
      
      for (unsigned int t = 1; t < obs.n_cols; t++) {
        alpha = transition.t() * alpha;
        for (unsigned int r = 0; r < obs.n_slices; r++) {
          alpha %= emission.slice(r).col(obs(k, t, r));
        }
        
        tmp = sum(alpha);
        ll(k) += log(tmp);
        alpha /= tmp;
      }
    }
    return ll;
}
