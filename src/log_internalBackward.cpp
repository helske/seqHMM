// internal backward algorithm using log-space

#include "forward_backward.h"
#include "logsumexp.h"
void log_internalBackward(const arma::mat& transition, const arma::cube& emission,
  const arma::ucube& obs, arma::cube& beta, unsigned int threads) {

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(beta, obs, emission,transition)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      beta.slice(k).col(obs.n_cols - 1).zeros();
      for (int t = (obs.n_cols - 2); t >= 0; t--) {
        arma::vec tmpbeta(transition.n_rows);
        for (unsigned int i = 0; i < transition.n_rows; i++) {
          tmpbeta = beta.slice(k).col(t + 1) + transition.row(i).t();
          for (unsigned int r = 0; r < obs.n_rows; r++) {
            tmpbeta += emission.slice(r).col(obs(r, t + 1, k));
          }
          beta(i, t, k) = logSumExp(tmpbeta);
        }
      }
    }
}
