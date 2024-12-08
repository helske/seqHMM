#include "state_obs_probs.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
Rcpp::List ame_obs_nhmm_singlechannel(
    arma::mat& eta_pi, arma::cube& eta_A, arma::cube& eta_B,
    const arma::umat& obs, const arma::uvec Ti,
    const bool icpt_only_pi1, const bool icpt_only_A1, const bool icpt_only_B1,
    const bool iv_A1, const bool iv_B1, const bool tv_A1, const bool tv_B1,
    const arma::mat& X1_pi, const arma::cube& X1_A,const arma::cube& X1_B,
    const bool icpt_only_pi2, const bool icpt_only_A2, const bool icpt_only_B2,
    const bool iv_A2, const bool iv_B2, const bool tv_A2, const bool tv_B2,
    const arma::mat& X2_pi, const arma::cube& X2_A,const arma::cube& X2_B,
    const arma::field<arma::mat>& boot_gamma_pi,
    const arma::field<arma::cube>& boot_gamma_A,
    const arma::field<arma::cube>& boot_gamma_B,
    const arma::uword start, const arma::vec& probs, const arma::umat& idx) {
  
  arma::uword S = eta_A.n_slices;
  arma::uword M = eta_B.n_rows + 1;
  arma::uword T = obs.n_rows;
  arma::uword N = obs.n_cols;
  nhmm_sc model1(
      S, X1_pi, X1_A, X1_B, Ti, icpt_only_pi1, icpt_only_A1,
      icpt_only_B1, iv_A1, iv_B1, tv_A1, tv_B1, obs, eta_pi, eta_A, eta_B
  );
  
  nhmm_sc model2(
      S, X2_pi, X2_A, X2_B, Ti, icpt_only_pi2, icpt_only_A2,
      icpt_only_B2, iv_A2, iv_B2, tv_A2, tv_B2, obs, eta_pi, eta_A, eta_B
  );
  arma::cube obs_prob1(M, T, N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob1(S, T, N, arma::fill::value(arma::datum::nan));
  model1.compute_state_obs_probs(start, obs_prob1, state_prob1);
  arma::cube obs_prob2(M, T, N, arma::fill::value(arma::datum::nan));
  arma::cube state_prob2(S, T, N, arma::fill::value(arma::datum::nan));
  model2.compute_state_obs_probs(start, obs_prob2, state_prob2);
  arma::mat point_estimate(M, T, arma::fill::value(arma::datum::nan));
  arma::mat diff(M, N);
  for (arma::uword t = start - 1; t < T; t++) {
    diff = obs_prob1.col(t) - obs_prob2.col(t);
    arma::uvec non_na = arma::find_finite(diff.row(0));
    if (!non_na.is_empty()) {
      point_estimate.col(t) = arma::mean(diff.cols(non_na), 1);
    }
  }
  arma::uword nsim = boot_gamma_pi.n_elem;
  if (nsim > 1) {
    arma::cube out(M, T, nsim, arma::fill::value(arma::datum::nan));
    for (arma::uword i = 0; i < nsim; i++) {
      model1.gamma_pi = boot_gamma_pi(i);
      for (arma::uword s = 0; s < S; s++) {
        model1.gamma_A.slice(s) = boot_gamma_A(i).slice(s);
        model1.gamma_B.slice(s) = boot_gamma_B(i).slice(s);
      }
      model1.compute_state_obs_probs(start, obs_prob1, state_prob1);
      model2.gamma_pi = model1.gamma_pi;
      model2.gamma_A = model1.gamma_A;
      model2.gamma_B = model1.gamma_B;
      model2.compute_state_obs_probs(start, obs_prob2, state_prob2);
      
      for (arma::uword t = start - 1; t < T; t++) {
        diff = obs_prob1.col(t) - obs_prob2.col(t);
        diff = diff.cols(idx.col(i));
        arma::uvec non_na = arma::find_finite(diff.row(0));
        if (!non_na.is_empty()) {
          out.slice(i).col(t) = arma::mean(diff.cols(non_na), 1);
        }
      }
    }
    
    arma::cube quantiles(M, T, probs.n_elem, arma::fill::value(arma::datum::nan));
    arma::mat tmp(M, nsim);
    for (arma::uword t = start - 1; t < T; t++) {
      tmp = out.col(t);
      quantiles.col(t) = arma::quantile(tmp, probs, 1);
    }
    return Rcpp::List::create(
      Rcpp::Named("point_estimate") = Rcpp::wrap(point_estimate),
      Rcpp::Named("quantiles") = Rcpp::wrap(quantiles)
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("point_estimate") = Rcpp::wrap(point_estimate)
  );
}

// [[Rcpp::export]]
Rcpp::List ame_obs_nhmm_multichannel(
    arma::mat& eta_pi, arma::cube& eta_A, arma::field<arma::cube>& eta_B,
    const arma::ucube& obs, const arma::uvec Ti,
    const bool icpt_only_pi1, const bool icpt_only_A1, const bool icpt_only_B1,
    const bool iv_A1, const bool iv_B1, const bool tv_A1, const bool tv_B1,
    const arma::mat& X1_pi, const arma::cube& X1_A, const arma::cube& X1_B,
    const bool icpt_only_pi2, const bool icpt_only_A2, const bool icpt_only_B2,
    const bool iv_A2, const bool iv_B2, const bool tv_A2, const bool tv_B2,
    const arma::mat& X2_pi, const arma::cube& X2_A,const arma::cube& X2_B,
    const arma::field<arma::mat>& boot_gamma_pi,
    const arma::field<arma::cube>& boot_gamma_A,
    const arma::field<arma::cube>& boot_gamma_B,
    const arma::uword start, const arma::vec& probs, const arma::umat& idx) {
  
  arma::uword S = eta_A.n_slices;
  arma::uword C = obs.n_rows;
  arma::uword T = obs.n_cols;
  arma::uword N = obs.n_slices;
  nhmm_mc model1(
      S, X1_pi, X1_A, X1_B, Ti, icpt_only_pi1, icpt_only_A1,
      icpt_only_B1, iv_A1, iv_B1, tv_A1, tv_B1, obs, eta_pi, eta_A, eta_B
  );
  
  nhmm_mc model2(
      S, X2_pi, X2_A, X2_B, Ti, icpt_only_pi2, icpt_only_A2,
      icpt_only_B2, iv_A2, iv_B2, tv_A2, tv_B2, obs, eta_pi, eta_A, eta_B
  );
  arma::uvec M = model1.M;
  
  arma::field<arma::cube> obs_prob1(C);
  for (arma::uword c = 0; c < C; c++) {
    obs_prob1(c) = arma::cube(M(c), T, N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob1(S, T, N, arma::fill::value(arma::datum::nan));
  model1.compute_state_obs_probs(start, obs_prob1, state_prob1);
  arma::field<arma::cube> obs_prob2(C);  
  for (arma::uword c = 0; c < C; c++) {
    obs_prob2(c) = arma::cube(M(c), T, N, arma::fill::value(arma::datum::nan));
  }
  arma::cube state_prob2(S, T, N, arma::fill::value(arma::datum::nan));
  model2.compute_state_obs_probs(start, obs_prob2, state_prob2);
  
  arma::field<arma::mat> point_estimate(C);
  arma::field<arma::mat> diff(C);
  for (arma::uword c = 0; c < C; c++) {
    point_estimate(c) = arma::mat(M(c), T, arma::fill::value(arma::datum::nan));
    diff(c) = arma::mat(M(c), N);
    
    for (arma::uword t = start - 1; t < T; t++) {
      diff(c) = obs_prob1(c).col(t) - obs_prob2(c).col(t);
      arma::uvec non_na = arma::find_finite(diff(c).row(0));
      if (!non_na.is_empty()) {
        point_estimate(c).col(t) = arma::mean(diff(c).cols(non_na), 1);
      }
    }
  }
  
  arma::uword nsim = boot_gamma_pi.n_elem;
  arma::field<arma::cube> out(C);
  for (arma::uword c = 0; c < C; c++) {
    out(c) = arma::cube(M(c), T, nsim, arma::fill::value(arma::datum::nan));
  }
  for (arma::uword i = 0; i < nsim; i++) {
    model1.gamma_pi = boot_gamma_pi(i);
    for (arma::uword s = 0; s < S; s++) {
      model1.gamma_A.slice(s) = boot_gamma_A(i).slice(s);
    }
    for (arma::uword c = 0; c < C; c++) {
      for (arma::uword s = 0; s < S; s++) {
        model1.gamma_B(c).slice(s) = boot_gamma_B(i, c).slice(s);
      }
    }
    model1.compute_state_obs_probs(start, obs_prob1, state_prob1);
    model2.gamma_pi = model1.gamma_pi;
    model2.gamma_A = model1.gamma_A;
    model2.gamma_B = model1.gamma_B;
    model2.compute_state_obs_probs(start, obs_prob2, state_prob2);
    for (arma::uword c = 0; c < C; c++) {
      for (arma::uword t = start - 1; t < T; t++) {
        diff(c) = obs_prob1(c).col(t) - obs_prob2(c).col(t);
        diff(c) = diff(c).cols(idx.col(i));
        arma::uvec non_na = arma::find_finite(diff(c).row(0));
        if (!non_na.is_empty()) {
          out(c).slice(i).col(t) = arma::mean(diff(c).cols(non_na), 1);
        }
      }
    }
  }
  
  arma::field<arma::cube> quantiles(C);
  for (arma::uword c = 0; c < C; c++) {
    quantiles(c) = arma::cube(M(c), T, probs.n_elem, arma::fill::value(arma::datum::nan));
    arma::mat tmp(M(c), nsim);
    for (arma::uword t = start - 1; t < T; t++) {
      tmp = out(c).col(t);
      quantiles(c).col(t) = arma::quantile(tmp, probs, 1);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("point_estimate") = Rcpp::wrap(point_estimate),
    Rcpp::Named("quantiles") = Rcpp::wrap(quantiles)
  );
}