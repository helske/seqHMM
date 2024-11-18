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

void nhmm_sc::compute_state_obs_probs(
    const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob) {
  
  if (start > 1) {
    obs_prob.cols(0, start - 2).fill(-arma::datum::inf);
  }
  bool not_updated = true;
  for (arma::uword i = 0; i < N; i++) {
    arma::uword upper_bound = std::min(start - 1, Ti(i));
    for (arma::uword t = 0; t < upper_bound; t++) {
      if (obs(t, i) < M) {
        obs_prob(obs(t, i), t, i) = 0;
      } else {
        obs_prob.slice(i).col(t).fill(arma::datum::nan); // could predict these as well
      }
    }
    if (start < Ti(i)) {
      if (!icpt_only_pi || not_updated) {
        update_pi(i);
      }
      if (iv_A || not_updated) {
        update_A(i);
      }
      if (iv_B || not_updated) {
        update_B(i);
      }
      not_updated = false;
      update_log_py(i);
      univariate_forward_nhmm(
        state_prob.slice(i), log_Pi, log_A, 
        log_py.cols(0, start - 1)
      );
      univariate_state_prob(
        state_prob.slice(i), start, Ti(i), log_A
      );
      univariate_obs_prob(
        obs_prob.slice(i), 
        state_prob.slice(i),
        start, Ti(i),
        log_B
      );
    }
  }
  obs_prob = arma::exp(obs_prob);
  state_prob = arma::exp(state_prob);
}

void nhmm_mc::compute_state_obs_probs(
    const arma::uword start, arma::field<arma::cube>& obs_prob, 
    arma::cube& state_prob) {
  
  if (start > 1) {
    for (arma::uword c = 0; c < C; c++) {
      obs_prob(c).cols(0, start - 2).fill(-arma::datum::inf);
    }
  }
  bool not_updated = true;
  for (arma::uword i = 0; i < N; i++) {
    arma::uword upper_bound = std::min(start - 1, Ti(i));
    for (arma::uword t = 0; t < upper_bound; t++) {
      for (arma::uword c = 0; c < C; c++) {
        if (obs(c, t, i) < M(c)) {
          obs_prob(c)(obs(c, t, i), t, i) = 0;
        } else {
          obs_prob(c).slice(i).col(t).fill(arma::datum::nan); // could predict these as well
        }
      }
    }
    if (start < Ti(i)) {
      if (!icpt_only_pi || not_updated) {
        update_pi(i);
      }
      if (iv_A || not_updated) {
        update_A(i);
      }
      if (iv_B || not_updated) {
        update_B(i);
      }
      not_updated = false;
      update_log_py(i);
      univariate_forward_nhmm(
        state_prob.slice(i), log_Pi, log_A, 
        log_py.cols(0, start - 1)
      );
      univariate_state_prob(
        state_prob.slice(i), start, Ti(i), log_A
      );
      for (arma::uword c = 0; c < C; c++) {
        univariate_obs_prob(
          obs_prob(c).slice(i), 
          state_prob.slice(i),
          start, Ti(i),
          log_B(c)
        );
      }
    }
  }
  for (arma::uword c = 0; c < C; c++) {
    obs_prob(c) = arma::exp(obs_prob(c));
  }
  state_prob = arma::exp(state_prob);
}

void mnhmm_sc::compute_state_obs_probs(
    const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob) {
  
  if (start > 1) {
    obs_prob.cols(0, start - 2).fill(-arma::datum::inf);
  }
  bool not_updated = true;
  arma::cube tmp(M, T, D);
  for (arma::uword i = 0; i < N; i++) {
    arma::uword upper_bound = std::min(start - 1, Ti(i));
    for (arma::uword t = 0; t < upper_bound; t++) {
      if (obs(t, i) < M) {
        obs_prob(obs(t, i), t, i) = 0;
      } else {
        obs_prob.slice(i).col(t).fill(arma::datum::nan); // could predict these as well
      }
    }
    if (start < Ti(i)) {
      if (!icpt_only_omega || not_updated) {
        update_omega(i);
      }
      if (!icpt_only_pi || not_updated) {
        update_pi(i);
      }
      if (iv_A || not_updated) {
        update_A(i);
      }
      if (iv_B || not_updated) {
        update_B(i);
      }
      not_updated = false;
      update_log_py(i);
      for (arma::uword d = 0; d < D; d++) {
        arma::subview<double> submat = 
          state_prob.slice(i).rows(d * S, (d + 1) * S - 1);
        univariate_forward_nhmm(
          submat,
          log_Pi(d),
          log_A(d), 
          log_py.slice(d).cols(0, start - 1)
        );
        univariate_state_prob(
          submat, start, Ti(i), log_A(d)
        );
        submat += log_omega(d);
        univariate_obs_prob(
          tmp.slice(d), 
          submat,
          start, Ti(i),
          log_B(d)
        );
      }
      for (arma::uword t = start - 1; t < Ti(i); t++) {
        for (arma::uword m = 0; m < M; m++) {
          obs_prob(m, t, i) = logSumExp(tmp.tube(m, t));
        }
      }
    }
  }
  obs_prob = arma::exp(obs_prob);
  state_prob = arma::exp(state_prob);
}

void mnhmm_mc::compute_state_obs_probs(
    const arma::uword start, arma::field<arma::cube>& obs_prob, 
    arma::cube& state_prob) {
  
  if (start > 1) {
    for (arma::uword c = 0; c < C; c++) {
      obs_prob(c).cols(0, start - 2).fill(-arma::datum::inf);
    }
  }
  bool not_updated = true;
  arma::field<arma::cube> tmp(C);
  for (arma::uword c = 0; c < C; c++) {
    tmp(c) = arma::cube(M(c), T, D);
  }
  
  for (arma::uword i = 0; i < N; i++) {
    arma::uword upper_bound = std::min(start - 1, Ti(i));
    for (arma::uword t = 0; t < upper_bound; t++) {
      for (arma::uword c = 0; c < C; c++) {
        if (obs(c, t, i) < M(c)) {
          obs_prob(c)(obs(c, t, i), t, i) = 0;
        } else {
          obs_prob(c).slice(i).col(t).fill(arma::datum::nan); // could predict these as well
        }
      }
    }
    
    if (start < Ti(i)) {
      if (!icpt_only_omega || not_updated) {
        update_omega(i);
      }
      if (!icpt_only_pi || not_updated) {
        update_pi(i);
      }
      if (iv_A || not_updated) {
        update_A(i);
      }
      if (iv_B || not_updated) {
        update_B(i);
      }
      not_updated = false;
      update_log_py(i);
      
      for (arma::uword d = 0; d < D; d++) {
        arma::subview<double> submat = 
          state_prob.slice(i).rows(d * S, (d + 1) * S - 1);
        univariate_forward_nhmm(
          submat,
          log_Pi(d),
          log_A(d), 
          log_py.slice(d).cols(0, start - 1)
        );
        
        univariate_state_prob(
          submat, start, Ti(i), log_A(d)
        );
        
        submat += log_omega(d);
        for (arma::uword c = 0; c < C; c++) {
          univariate_obs_prob(
            tmp(c).slice(d), 
            submat,
            start, Ti(i),
            log_B(c, d)
          );
        }
      }
      for (arma::uword c = 0; c < C; c++) {
        for (arma::uword t = start - 1; t < Ti(i); t++) {
          for (arma::uword m = 0; m < M(c); m++) {
            obs_prob(c)(m, t, i) = logSumExp(tmp(c).tube(m, t));
          }
        }
      }
    }
  }
  for (arma::uword c = 0; c < C; c++) {
    obs_prob(c) = arma::exp(obs_prob(c));
  }
  state_prob = arma::exp(state_prob);
}
#endif

