#include "seqHMM.h"
// [[Rcpp::export]]

NumericVector logLikHMM(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);

  NumericVector ll(obs.n_rows);

#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition)
  for (int k = 0; k < obs.n_rows; k++) {
    arma::vec alpha = init;
    
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

