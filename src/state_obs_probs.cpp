// state_obs_probs algorithm for NHMM
#include "state_obs_probs.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
Rcpp::List state_obs_probs_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_nhmm_multichannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti,
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B,
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A,
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    obs_prob(c) = arma::cube(model.M(c), model.T, model.N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_mnhmm_singlechannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti,
    const bool icpt_only_omega, const bool icpt_only_pi,
    const bool icpt_only_A, const bool icpt_only_B, const bool iv_A,
    const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  mnhmm_sc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti,
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob(model.S * model.D, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
}

// [[Rcpp::export]]
Rcpp::List state_obs_probs_mnhmm_multichannel(
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const arma::uvec Ti,
    const bool icpt_only_omega, const bool icpt_only_pi,
    const bool icpt_only_A, const bool icpt_only_B, const bool iv_A,
    const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start) {
  
  mnhmm_mc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti,
      icpt_only_omega, icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    obs_prob(c) = arma::cube(model.M(c), model.T, model.N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob(model.S * model.D, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs(start, obs_prob, state_prob);
  
  return Rcpp::List::create(
    Rcpp::Named("obs_prob") = Rcpp::wrap(obs_prob), 
    Rcpp::Named("state_prob") = Rcpp::wrap(state_prob)
  );
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
    if (start <= Ti(i)) {
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
      univariate_state_prob(
        state_prob.slice(i), start, Ti(i), log_pi, log_A, log_py
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
    if (start <= Ti(i)) {
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
      univariate_state_prob(
        state_prob.slice(i), start, Ti(i), log_pi, log_A, log_py
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
    if (start <= Ti(i)) {
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
        univariate_state_prob(
          submat, start, Ti(i), log_pi(d), log_A(d), log_py.slice(d)
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
    
    if (start <= Ti(i)) {
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
        univariate_state_prob(
          submat, start, Ti(i), log_pi(d), log_A(d), log_py.slice(d)
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
