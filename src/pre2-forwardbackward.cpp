#include "pre2-forward_backward.h"
#include "reparma.h"
#include "useomp.h"

// Forward-backward algorithm for non-mixture hidden Markov models
// [[Rcpp::export]]
Rcpp::List forwardbackward(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, bool forwardonly, arma::uword threads) {
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //n,k
  
  internalForward(transition.t(), emission, init, obs, alpha, scales, threads);
  
  if(!scales.is_finite()) {
    Rcpp::stop("Scaling factors contain non-finite values. \n Check the model or try using the log-space version of the algorithm.");
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }
  
  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha),
      Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    internalBackward(transition, emission, obs, beta, scales, threads);
    if(!beta.is_finite()) {
      Rcpp::stop("Backward probabilities contain non-finite values. Check the model or try using the log-space version of the algorithm.");
    }
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta),
      Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  }
}

// Forward-backward algorithm for mixture hidden Markov models
// [[Rcpp::export]]
Rcpp::List forwardbackwardx(const arma::mat& transition, const arma::cube& emission,
                            const arma::vec& init, const arma::ucube obs, const arma::mat& coef, const arma::mat& X,
                            const arma::uvec& numberOfStates, bool forwardonly, arma::uword threads) {
  
  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= arma::sum(weights, 0);
  
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //m,n,k
  
  internalForward(transition.t(), emission, initk, obs, alpha, scales, threads);
  
  if(!scales.is_finite()) {
    Rcpp::stop("Scaling factors contain non-finite values. Check the model or try using the log-space version of the algorithm.");
  }
  double max_sf = scales.max();
  if (max_sf > 1e150) {
    Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
  }
  
  if (forwardonly) {
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha),
                              Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  } else {
    arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
    internalBackward(transition, emission, obs, beta, scales, threads);
    if(!beta.is_finite()) {
      Rcpp::stop("Backward probabilities contain non-finite values. Check the model or try using the log-space version of the algorithm.");
    }
    return Rcpp::List::create(Rcpp::Named("forward_probs") = Rcpp::wrap(alpha), Rcpp::Named("backward_probs") = Rcpp::wrap(beta),
                              Rcpp::Named("scaling_factors") = Rcpp::wrap(scales));
  }
  
  return Rcpp::wrap(alpha);
}

void internalForward(const arma::mat& transition_t, const arma::cube& emission, const arma::vec& init,
                     const arma::ucube& obs, arma::cube& alpha, arma::mat& scales, arma::uword threads) {
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(alpha, scales, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      alpha.slice(k).col(0) = init;
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; r++) {
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
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      
      alpha.slice(k).col(0) = init.col(k);
      for (arma::uword r = 0; r < obs.n_rows; r++) {
        alpha.slice(k).col(0) %= emission.slice(r).col(obs(r, 0, k));
      }
      scales(0, k) = 1.0 / sum(alpha.slice(k).col(0));
      alpha.slice(k).col(0) *= scales(0, k);
      for (arma::uword t = 1; t < obs.n_cols; t++) {
        alpha.slice(k).col(t) = transition_t * alpha.slice(k).col(t - 1);
        for (arma::uword r = 0; r < obs.n_rows; r++) {
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
    for (arma::uword k = 0; k < obs.n_slices; k++) {
      beta.slice(k).col(obs.n_cols - 1).fill(scales(obs.n_cols - 1, k));
      for (int t = obs.n_cols - 2; t >= 0; t--) {
        arma::vec tmpbeta = beta.slice(k).col(t + 1);
        for (arma::uword r = 0; r < obs.n_rows; r++) {
          tmpbeta %= emission.slice(r).col(obs(r, t + 1, k));
        }
        beta.slice(k).col(t) = transition * tmpbeta * scales(t, k);
      }
    }
}


void uvForward(const arma::mat& transition_t, const arma::cube& emission, const arma::vec& init,
               const arma::umat& obs, arma::mat& alpha, arma::vec& scales) {
  
  alpha.col(0) = init;
  for (arma::uword r = 0; r < obs.n_rows; r++) {
    alpha.col(0) %= emission.slice(r).col(obs(r, 0));
  }
  scales(0) = 1.0 / sum(alpha.col(0));
  alpha.col(0) *= scales(0);
  for (arma::uword t = 1; t < obs.n_cols; t++) {
    alpha.col(t) = transition_t * alpha.col(t - 1);
    for (arma::uword r = 0; r < obs.n_rows; r++) {
      alpha.col(t) %= emission.slice(r).col(obs(r, t));
    }
    scales(t) = 1.0 / sum(alpha.col(t));
    alpha.col(t) *= scales(t);
  }
  
}

void uvBackward(const arma::mat& transition, const arma::cube& emission,
                const arma::umat& obs, arma::mat& beta, const arma::vec& scales) {
  
  beta.col(obs.n_cols - 1).fill(scales(obs.n_cols - 1));
  for (int t = obs.n_cols - 2; t >= 0; t--) {
    arma::vec tmpbeta = beta.col(t + 1);
    for (arma::uword r = 0; r < obs.n_rows; r++) {
      tmpbeta %= emission.slice(r).col(obs(r, t + 1));
    }
    beta.col(t) =  transition * tmpbeta * scales(t);
  }
}
