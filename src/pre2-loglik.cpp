#include "logsumexp.h"
#include "reparma.h"
#include "useomp.h"
// log-likelihood of HMM
// [[Rcpp::export]]
Rcpp::NumericVector logLikHMM(const arma::mat& transition, const arma::cube& emission,
                              const arma::vec& init, const arma::ucube& obs, arma::uword threads) {
  
  arma::vec ll(obs.n_slices);
  arma::mat transition_t(transition.t());
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition_t)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      arma::vec alpha = init;
      
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha %= emission.slice(r).col(obs(r, 0, k));
      }
      
      double tmp = sum(alpha);
      ll(k) = log(tmp);
      alpha /= tmp;
      
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        alpha = transition_t * alpha;
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha %= emission.slice(r).col(obs(r, t, k));
        }
        
        tmp = sum(alpha);
        ll(k) += log(tmp);
        alpha /= tmp;
      }
    }
    return Rcpp::wrap(ll);
}


// log-likelihood of MHMM
// [[Rcpp::export]]
Rcpp::NumericVector logLikMixHMM(const arma::mat& transition, const arma::cube& emission,
                                 const arma::vec& init, const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
                                 const arma::uvec& numberOfStates, arma::uword threads) {
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::wrap(-arma::datum::inf);
  }
  weights.each_row() /= sum(weights, 0);
  
  arma::vec ll(obs.n_slices);
  arma::mat transition_t(transition.t());
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition_t, numberOfStates)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      arma::vec alpha = init % reparma(weights.col(k), numberOfStates);
      
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha %= emission.slice(r).col(obs(r, 0, k));
      }
      
      double tmp = sum(alpha);
      ll(k) = log(tmp);
      alpha /= tmp;
      
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        alpha = transition_t * alpha;
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          alpha %= emission.slice(r).col(obs(r, t, k));
        }
        
        tmp = sum(alpha);
        ll(k) += log(tmp);
        alpha /= tmp;
      }
    }
    return Rcpp::wrap(ll);
}

// log-likelihood of HMM using log-space
// [[Rcpp::export]]
Rcpp::NumericVector log_logLikHMM(const arma::mat& transition_, const arma::cube& emission_, 
                                  const arma::vec& init_, const arma::ucube& obs, arma::uword threads) {
  
  arma::vec init = log(init_);
  arma::mat transition = log(transition_);
  arma::cube emission = log(emission_);
  
  arma::vec ll(obs.n_slices);
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, init, emission, transition)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      arma::vec alpha = init;
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha += emission.slice(r).col(obs(r, 0, k));
      }
      
      arma::vec alphatmp(emission.n_rows);
      
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        for (arma::uword i = 0; i < emission.n_rows; ++i) {
          alphatmp(i) = logSumExp(alpha + transition.col(i));
          for (arma::uword r = 0; r < obs.n_rows; ++r) {
            alphatmp(i) += emission(i, obs(r, t, k), r);
          }
        }
        alpha = alphatmp;
      }
      ll(k) = logSumExp(alpha);
    }
    return Rcpp::wrap(ll);
}

// log-likelihood of MHMM using log-space
// [[Rcpp::export]]
Rcpp::NumericVector log_logLikMixHMM(arma::mat transition, arma::cube emission, arma::vec init,
                                     const arma::ucube& obs, const arma::mat& coef, const arma::mat& X,
                                     const arma::uvec& numberOfStates, arma::uword threads) {
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::wrap(-arma::datum::inf);
  }
  weights.each_row() /= sum(weights, 0);
  
  weights = log(weights);
  transition = log(transition);
  emission = log(emission);
  init = log(init);
  
  arma::vec ll(obs.n_slices);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(ll, obs, weights, init, emission, transition, numberOfStates)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      arma::vec alpha = init + reparma(weights.col(k), numberOfStates);
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        alpha += emission.slice(r).col(obs(r, 0, k));
      }
      arma::vec alphatmp(emission.n_rows);
      for (arma::uword t = 1; t < obs.n_cols; ++t) {
        for (arma::uword i = 0; i < emission.n_rows; ++i) {
          alphatmp(i) = logSumExp(alpha + transition.col(i));
          for (arma::uword r = 0; r < obs.n_rows; ++r) {
            alphatmp(i) += emission(i, obs(r, t, k), r);
          }
        }
        alpha = alphatmp;
      }
      ll(k) = logSumExp(alpha);
    }
    return Rcpp::wrap(ll);
}


