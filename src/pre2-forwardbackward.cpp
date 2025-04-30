#include "pre2-forward_backward.h"
#include "reparma.h"
#include "useomp.h"

void internalForward(const arma::mat& transition_t, const arma::cube& emission, const arma::vec& init,
                     const arma::ucube& obs, arma::cube& alpha, arma::mat& scales, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      alpha.slice(k).col(0) = init;
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha.slice(k).col(t) %= emission.slice(r).col(obs(r, t, k));
        }
        scales(t, k) = 1.0 / sum(alpha.slice(k).col(t));
        alpha.slice(k).col(t) *= scales(t, k);
      }
    }
}

void internalForward(const arma::mat& transition_t, const arma::cube& emission,
                      const arma::mat& init, const arma::ucube& obs, arma::cube& alpha, arma::mat& scales,
                      arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      
      alpha.slice(k).col(0) = init.col(k);
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha.slice(k).col(t) %= emission.slice(r).col(obs(r, t, k));
        }
        scales(t, k) = 1.0 / sum(alpha.slice(k).col(t));
        alpha.slice(k).col(t) *= scales(t, k);
      }
    }
}
void internalBackward(const arma::mat& transition, const arma::cube& emission,
                      const arma::ucube& obs, arma::cube& beta, const arma::mat& scales, 
                      arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(beta, scales, obs, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      beta.slice(k).col(obs.n_cols - 1).fill(scales(obs.n_cols - 1, k));
      for (arma::uword t = obs.n_cols - 1; t-- > 0;) {
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          tmpbeta %= emission.slice(r).col(obs(r, t + 1, k));
        }
        beta.slice(k).col(t) = transition * tmpbeta * scales(t, k);
      }
    }
}


void uvForward(const arma::mat& transition_t, const arma::cube& emission, const arma::vec& init,
               const arma::umat& obs, arma::mat& alpha, arma::vec& scales) {
  
  alpha.col(0) = init;
  for (arma::uword r = 0; r < obs.n_rows; ++r) {
    alpha.col(0) %= emission.slice(r).col(obs(r, 0));
  }
  scales(0) = 1.0 / sum(alpha.col(0));
  alpha.col(0) *= scales(0);
  for (arma::uword t = 1; t < obs.n_cols; ++t) {
    alpha.col(t) = transition_t * alpha.col(t - 1);
    for (arma::uword r = 0; r < obs.n_rows; ++r) {
      alpha.col(t) %= emission.slice(r).col(obs(r, t));
    }
    scales(t) = 1.0 / sum(alpha.col(t));
    alpha.col(t) *= scales(t);
  }
  
}

void uvBackward(const arma::mat& transition, const arma::cube& emission,
                const arma::umat& obs, arma::mat& beta, const arma::vec& scales) {
  
  beta.col(obs.n_cols - 1).fill(scales(obs.n_cols - 1));
  for (arma::uword t = obs.n_cols - 1; t-- > 0;) {
    arma::vec tmpbeta = beta.col(t + 1);
    for (arma::uword r = 0; r < obs.n_rows; ++r) {
      tmpbeta %= emission.slice(r).col(obs(r, t + 1));
    }
    beta.col(t) =  transition * tmpbeta * scales(t);
  }
}
