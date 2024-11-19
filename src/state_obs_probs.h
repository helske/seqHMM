#ifndef AME_OBS_NHMM_H
#define AME_OBS_NHMM_H

#include <RcppArmadillo.h>
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"
#include "nhmm_forward.h"
#include "logsumexp.h"

template<typename submat>
void univariate_state_prob(
    submat& log_state_prob,
    const arma::uword start, 
    const arma::uword end, 
    const arma::cube& log_A) {
  
  arma::uword S = log_A.n_rows;
  for (arma::uword t = 0; t < start; t++) {
    log_state_prob.col(t) -= logSumExp(log_state_prob.col(t));
  }
  for (arma::uword t = start; t < end; t++) {
    for (arma::uword s = 0; s < S; s++) {
      log_state_prob(s, t) = logSumExp(log_state_prob.col(t - 1) + log_A.slice(t - 1).col(s));
    }
    log_state_prob.col(t) -= logSumExp(log_state_prob.col(t));
  }
}

template<typename submat>
void univariate_obs_prob(
    arma::mat& log_obs_prob,
    const submat& log_state_prob,
    const arma::uword start,
    const arma::uword end, 
    const arma::cube& log_B) {
  arma::uword M = log_B.n_cols - 1;
  for (arma::uword t = start - 1; t < end; t++) {
    for (arma::uword m = 0; m < M; m++) {
      log_obs_prob(m, t) = logSumExp(log_state_prob.col(t) + 
        log_B.slice(t).col(m));
    }
  }
}
#endif

