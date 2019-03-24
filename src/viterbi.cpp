//Viterbi algorithm for HMM
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List viterbi(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs) {

  arma::umat q(obs.n_slices, obs.n_cols);
  arma::vec logp(obs.n_slices);

  arma::mat delta(emission.n_rows, obs.n_cols);
  arma::umat phi(emission.n_rows, obs.n_cols);

  for (unsigned int k = 0; k < obs.n_slices; k++) {

    delta.col(0) = init;
    for (unsigned int r = 0; r < emission.n_slices; r++) {
      delta.col(0) += emission.slice(r).col(obs(r, 0, k));
    }

    phi.col(0).zeros();

    for (unsigned int t = 1; t < obs.n_cols; t++) {
      for (unsigned int j = 0; j < emission.n_rows; j++) {
        phi(j, t) = (delta.col(t - 1) + transition.col(j)).index_max();
        delta(j, t) = delta(phi(j, t), t - 1) + transition(phi(j, t), j);
        for (unsigned int r = 0; r < emission.n_slices; r++) {
          delta(j, t) += emission(j, obs(r, t, k), r);
        }
      }
    }
    
   // Changes in Armadillo
   // delta.col(obs.n_cols - 1).max(q(k, obs.n_cols - 1));
   q(k, obs.n_cols - 1) = delta.col(obs.n_cols - 1).index_max();
    for (int t = (obs.n_cols - 2); t >= 0; t--) {
      q(k, t) = phi(q(k, t + 1), t + 1);
    }
    logp(k) = delta.col(obs.n_cols - 1).max();
  }

  return Rcpp::List::create(Rcpp::Named("q") = Rcpp::wrap(q), Rcpp::Named("logp") = Rcpp::wrap(logp));
}
