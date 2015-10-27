#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector log_logLikHMM(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, int threads = 1) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);

  transition = log(transition);
  emission = log(emission);
  init = log(init);

  NumericVector ll(obs.n_rows);
#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition)
  for (int k = 0; k < obs.n_rows; k++) {
    arma::vec alpha(emission.n_rows);
    for (int i = 0; i < emission.n_rows; i++) {
      alpha(i) = init(i);
      for (int r = 0; r < obs.n_slices; r++) {
        alpha(i) += emission(i, obs(k, 0, r), r);
      }
    }

    arma::vec alphatmp(emission.n_rows);

    for (int t = 1; t < obs.n_cols; t++) {
      for (int i = 0; i < emission.n_rows; i++) {
        alphatmp(i) = logSumExp(alpha + transition.col(i));
        for (int r = 0; r < obs.n_slices; r++) {
          alphatmp(i) += emission(i, obs(k, t, r), r);
        }
      }
      alpha = alphatmp;
    }
    ll(k) = logSumExp(alpha);
  }

  return ll;
}

