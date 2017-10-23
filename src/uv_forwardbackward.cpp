// univariate forward and backward algorithms for HMMs and MHMMs

#include "forward_backward.h"

void uvForward(const arma::sp_mat& transition_t, const arma::cube& emission, const arma::vec& init,
  const arma::umat& obs, arma::mat& alpha, arma::vec& scales) {

  alpha.col(0) = init;
  for (unsigned int r = 0; r < obs.n_rows; r++) {
    alpha.col(0) %= emission.slice(r).col(obs(r, 0));
  }
  scales(0) = 1.0 / sum(alpha.col(0));
  alpha.col(0) *= scales(0);
  for (unsigned int t = 1; t < obs.n_cols; t++) {
    alpha.col(t) = transition_t * alpha.col(t - 1);
    for (unsigned int r = 0; r < obs.n_rows; r++) {
      alpha.col(t) %= emission.slice(r).col(obs(r, t));
    }
    scales(t) = 1.0 / sum(alpha.col(t));
    alpha.col(t) *= scales(t);
  }

}

void uvBackward(const arma::sp_mat& transition, const arma::cube& emission,
  const arma::umat& obs, arma::mat& beta, const arma::vec& scales) {

  beta.col(obs.n_cols - 1).fill(scales(obs.n_cols - 1));
  for (int t = obs.n_cols - 2; t >= 0; t--) {
    arma::vec tmpbeta = beta.col(t + 1);
    for (unsigned int r = 0; r < obs.n_rows; r++) {
      tmpbeta %= emission.slice(r).col(obs(r, t + 1));
    }
    beta.col(t) =  transition * tmpbeta * scales(t);
  }
}

