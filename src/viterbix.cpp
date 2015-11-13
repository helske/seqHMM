#include "seqHMM.h"

// [[Rcpp::export]]

List viterbix(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, NumericMatrix coefs, NumericMatrix X_,
    IntegerVector numberOfStates) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);

  arma::umat q(obs.n_rows, obs.n_cols);
  arma::vec logp(obs.n_rows);


  int qn = coefs.nrow();
  arma::mat coef(coefs.begin(), qn, numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(), obs.n_rows, qn);

  arma::mat lweights = exp(X * coef).t();
  lweights.each_row() /= sum(lweights, 0);
  lweights = log(lweights);

  for (unsigned int k = 0; k < obs.n_rows; k++) {
    arma::mat delta(emission.n_rows, obs.n_cols);
    arma::umat phi(emission.n_rows, obs.n_cols);
    
    delta.col(0) = init + reparma(lweights.col(k), numberOfStates);
    for (unsigned int r = 0; r < emission.n_slices; r++) {
      delta.col(0) += emission.slice(r).col(obs(k, 0, r));
    }

    phi.col(0).zeros();

    for (unsigned int t = 1; t < obs.n_cols; t++) {
      for (unsigned int j = 0; j < emission.n_rows; j++) {
        (delta.col(t - 1) + transition.col(j)).max(phi(j, t));
        delta(j, t) = delta(phi(j, t), t - 1) + transition(phi(j, t), j);
        for (unsigned int r = 0; r < emission.n_slices; r++) {
          delta(j, t) += emission(j, obs(k, t, r), r);
        }
      }
    }

    delta.col(obs.n_cols - 1).max(q(k, obs.n_cols - 1));

    for (int t = (obs.n_cols - 2); t >= 0; t--) {
      q(k, t) = phi(q(k, t + 1), t + 1);
    }
    logp(k) = delta.col(obs.n_cols - 1).max();
  }

  return List::create(Named("q") = wrap(q), Named("logp") = wrap(logp));
}
