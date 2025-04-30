
#include "pre2-log_forward_backward.h"
#include "reparma.h"
#include "logsumexp.h"
#include "useomp.h"
// Forward-backward algorithm for non-mixture hidden Markov models using log-space
// [[Rcpp::export]]
Rcpp::List log_forwardbackward(const arma::mat& transition_, const arma::cube& emission_, 
  const arma::vec& init_, const arma::ucube& obs, bool forwardonly, arma::uword threads) {

  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k

  log_internalForward(transition, emission, init, obs, alpha, threads);

  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta));
  }

}

// [[Rcpp::export]]
Rcpp::List log_forwardbackwardx(const arma::mat& transition_, 
  const arma::cube& emission_, const arma::vec& init_,
  const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
    const arma::uvec& numberOfStates, bool forwardonly, arma::uword threads) {

  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);
  
  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= arma::sum(weights, 0);
  weights = log(weights);

  arma::mat initk(emission.n_rows, obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    initk.col(k) = init + reparma(weights.col(k), numberOfStates);
  }
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  log_internalForward(transition, emission, initk, obs, alpha, threads);

  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    log_internalBackward(transition, emission, obs, beta, threads);
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), 
                              Rcpp::Named("backward_probs") = Rcpp::wrap(beta));
  }

  return Rcpp::wrap(alpha);
}

// Internal forward algorithm for HMMs using log-space
void log_internalForward(const arma::mat& transition, const arma::cube& emission,
                         const arma::vec& init, const arma::ucube& obs, 
                         arma::cube& alpha, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      alpha.slice(k).col(0) = init;
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        for (arma::uword i = 0; i < transition.n_rows; ++i) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}

void log_internalForward(const arma::mat& transition, const arma::cube& emission,
                         const arma::mat& init, const arma::ucube& obs, 
                         arma::cube& alpha, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, obs, init, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      alpha.slice(k).col(0) = init.col(k);
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha.slice(k).col(0) += emission.slice(r).col(obs(r, 0, k));
      }
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        for (arma::uword i = 0; i < transition.n_rows; ++i) {
          alpha(i, t, k) = logSumExp(alpha.slice(k).col(t - 1) + transition.col(i));
        }
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha.slice(k).col(t) += emission.slice(r).col(obs(r, t, k));
        }
      }
    }
}
// internal backward algorithm using log-space
void log_internalBackward(const arma::mat& transition, const arma::cube& emission,
                          const arma::ucube& obs, arma::cube& beta, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(beta, obs, emission,transition)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      beta.slice(k).col(obs.n_cols - 1).zeros();
      for (arma::uword t = obs.n_cols - 1; t-- > 0;) {
        arma::vec tmpbeta(transition.n_rows);
        for (arma::uword i = 0; i < transition.n_rows; ++i) {
          tmpbeta = beta.slice(k).col(t + 1) + transition.row(i).t();
          for (arma::uword r = 0; r < obs.n_rows; ++r) {
            tmpbeta += emission.slice(r).col(obs(r, t + 1, k));
          }
          beta(i, t, k) = logSumExp(tmpbeta);
        }
      }
    }
}
