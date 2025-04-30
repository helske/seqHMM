//Viterbi algorithm for NHMM and MHMM, single_sequence
#include "viterbi.h"

double univariate_viterbi(
    arma::uvec& q,
    const arma::vec& log_pi, 
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  arma::mat delta(S, T);
  arma::umat phi(S, T);
  delta.col(0) = log_pi + log_py.col(0);
  phi.col(0).zeros();
  for (arma::uword t = 1; t < T; ++t) {
    for (arma::uword j = 0; j < S; ++j) {
      phi(j, t) = (delta.col(t - 1) + log_A.slice(t).col(j)).index_max();
      delta(j, t) = delta(phi(j, t), t - 1) + log_A(phi(j, t), j, t) + log_py(j, t);
    }
  }
  q(T - 1) = delta.col(T - 1).index_max();
  for (arma::uword t = (T - 1); t-- > 0;) {
    q(t) = phi(q(t + 1), t + 1);
  }
  return delta.col(T - 1).max();
}
