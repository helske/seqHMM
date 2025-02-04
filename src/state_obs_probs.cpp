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


// [[Rcpp::export]]
Rcpp::List state_obs_probs_fanhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uword start,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob(model.S, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.compute_state_obs_probs_fanhmm(start, obs_prob, state_prob, W_A, W_B);
  
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
    arma::uword upper_bound = std::min(start, Ti(i)) - 1;
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
    arma::uword upper_bound = std::min(start, Ti(i)) - 1;
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
    arma::uword upper_bound = std::min(start, Ti(i)) - 1;
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
    arma::uword upper_bound = std::min(start, Ti(i)) - 1;
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

void nhmm_sc::compute_state_obs_probs_fanhmm(
    const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B) {
  
  if (start > 1) {
    obs_prob.cols(0, start - 2).fill(-arma::datum::inf);
  }
  bool not_updated = true;
  arma::cube log_A_tm1(S, S, M);
  arma::cube log_B_t(S, M, M);
  for (arma::uword i = 0; i < N; i++) {
    arma::uword upper_bound = std::min(start, Ti(i)) - 1;
    for (arma::uword t = 0; t < upper_bound; t++) {
      obs_prob(obs(t, i), t, i) = 0;
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
     
      state_prob.slice(i).col(0) = log_pi;
      // forward algorithm until start - 1
      if (start > 1) {
        state_prob.slice(i).col(0) += log_py.col(0);
        for (arma::uword t = 1; t < start - 1; t++) {
          for (arma::uword s = 0; s < S; s++) {
            state_prob(s, t, i) = logSumExp(
              state_prob.slice(i).col(t - 1) + log_A.slice(t - 1).col(s) + log_py(s, t)
            );
          }
        }
        
        // predict one step ahead (start)
        for (arma::uword s = 0; s < S; s++) {
          state_prob(s, start - 1, i) = logSumExp(
            state_prob.slice(i).col(start - 2) + log_A.slice(start - 2).col(s)
          );
        }
      
        // normalize all
        for (arma::uword t = 0; t < start; t++) {
          state_prob.slice(i).col(t) -= logSumExp(state_prob.slice(i).col(t));
        }
        
       
        // y_start, no need to marginalize over missing as y_start-1 is still observed
        for (arma::uword m = 0; m < M; m++) {
          obs_prob(m, start - 1, i) = logSumExp(state_prob.slice(i).col(start - 1) +
            log_B.slice(start - 1).col(m));
        }
      }
      // // predict t = start, start-1 is the last observation (used in A_t-1 and B_t)
      // for (arma::uword s = 0; s < S; s++) {
      //   state_prob(s, start, i) = logSumExp(
      //     state_prob.slice(i).col(start - 1) + log_A.slice(start - 1).col(s)
      //   );
      // }
      // state_prob.slice(i).col(start) -= logSumExp(state_prob.slice(i).col(start));
      // for (arma::uword m = 0; m < M; m++) {
      //   obs_prob(m, start, i) = logSumExp(state_prob.slice(i).col(start) +
      //     log_B.slice(start).col(m));
      // }
      // predict t = start + 1 by marginalizing missing y_t-1
      arma::mat log_state_prob_new(S, M);
      arma::mat log_state_prob_old(S, M, arma::fill::value(-arma::datum::inf));
      log_state_prob_old.col(obs(start - 1, i)) = state_prob.slice(i).col(start - 1);
      state_prob.slice(i).cols(0, start - 1) = arma::exp(state_prob.slice(i).cols(0, start - 1));
      obs_prob.slice(i).cols(0, start - 1) = arma::exp(obs_prob.slice(i).cols(0, start - 1));
      arma::mat tmp(S, M);
      for (arma::uword t = start; t < Ti(i); t++) {
        // transition and emission matrices
        for (arma::uword m = 0; m < M; m++) {
          for (arma::uword s = 0; s < S; s ++) { // from states
            log_A_tm1.slice(m).row(s) = arma::log(softmax(gamma_A.slice(s) * W_A(m).slice(i).col(t - start))).t();
            log_B_t.slice(m).row(s) = arma::log(softmax(gamma_B.slice(s) * W_B(m).slice(i).col(t - start + 1))).t();
          }
        }
        for (arma::uword m = 0; m < M; m++) {
          for (arma::uword s = 0; s < S; s++) {
            for (arma::uword j = 0; j < M; j++) {
              for (arma::uword r = 0; r < S; r++) {
                tmp(r, j) = log_state_prob_old(r, j) + log_A_tm1(r, s, j) + log_B_t(s, m, j);
              }
            }
            // forward variable alpha(s, m) in log-space
            log_state_prob_new(s, m) = logSumExp(arma::vectorise(tmp));
          }
        }
        log_state_prob_old = log_state_prob_new;
        log_state_prob_new = arma::exp(log_state_prob_new);
        state_prob.slice(i).col(t) = arma::sum(log_state_prob_new, 1);
        obs_prob.slice(i).col(t) = arma::sum(log_state_prob_new, 0).t();
      }
    }
  }
}
