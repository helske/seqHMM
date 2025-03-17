// forward predictions for NHMMs
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
arma::field<arma::cube>  predict_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; i++) {
    obs_prob(i) = arma::cube(model.S, model.M, model.T, arma::fill::value(arma::datum::nan));
  }
  model.predict(obs_prob);
  return obs_prob;
}

// [[Rcpp::export]]
arma::field<arma::cube> predict_fanhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const bool autoregression) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; i++) {
    obs_prob(i) = arma::cube(model.S, model.M, model.T, arma::fill::value(arma::datum::nan));
  }
  model.predict_fanhmm(obs_prob, W_A, W_B, autoregression);
  
  return obs_prob;
}
// [[Rcpp::export]]
arma::field<arma::cube>  boot_predict_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::mat& gamma_pi, const arma::cube& gamma_A,
    const arma::cube& gamma_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; i++) {
    obs_prob(i) = arma::cube(model.S, model.M, model.T, arma::fill::value(arma::datum::nan));
  }
  model.predict(obs_prob);
  return obs_prob;
}

// [[Rcpp::export]]
arma::field<arma::cube> boot_predict_fanhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const arma::mat& gamma_pi, const arma::cube& gamma_A,
    const arma::cube& gamma_B, const bool autoregression) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; i++) {
    obs_prob(i) = arma::cube(model.S, model.M, model.T, arma::fill::value(arma::datum::nan));
  }
  model.predict_fanhmm(obs_prob, W_A, W_B, autoregression);
  return obs_prob;
}

void nhmm_sc::predict(arma::field<arma::cube>& obs_prob) {
  
  arma::vec alpha(S);
  for (arma::uword i = 0; i < N; i++) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (iv_B || i == 0) {
      update_B(i);
    }
    alpha = pi;
    obs_prob(i).slice(0) = B.slice(0).cols(0, M - 1);
    obs_prob(i).slice(0).each_col() %= alpha;
    alpha = obs_prob(i).slice(0).col(obs(0, i));
    alpha /= arma::accu(alpha);
    for (arma::uword t = 1; t < Ti(i); t++) {
      alpha = A.slice(t).t() * alpha;
      obs_prob(i).slice(t) = B.slice(t).cols(0, M - 1);
      obs_prob(i).slice(t).each_col() %= alpha;
      alpha = obs_prob(i).slice(t).col(obs(t, i));
      alpha /= arma::accu(alpha);
    }
  }
}

void nhmm_sc::predict_fanhmm(
    arma::field<arma::cube>& obs_prob, 
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B, 
    const bool autoregression) {
  
  arma::mat A_t(S, S);
  arma::mat B_t(M, S);
  arma::mat alpha_new(S, M); // P(alpha_t = s, y_t = m)
  arma::mat Btmp(S, M);
  arma::vec alpha(S);
  for (arma::uword i = 0; i < N; i++) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (iv_B || i == 0) {
      update_B(i);
    }
    // P(z_1)
    alpha = pi;
    // P(y_1)
    if (autoregression) {
      // first observation is fixed
      obs_prob(i).slice(0).zeros();
      obs_prob(i).slice(0).col(obs(0, i)) = alpha;
    } else {
      // no autoregressive component
      obs_prob(i).slice(0) = B.slice(0).cols(0, M - 1);
      obs_prob(i).slice(0).each_col() %= alpha;
      if (obs(0, i) < M) {
        alpha = obs_prob(i).slice(0).col(obs(0, i));
        alpha /= arma::accu(alpha);
      }
    }
    for (arma::uword t = 1; t < Ti(i); t++) {
      if (obs(t - 1, i) < M) {
        alpha = A.slice(t).t() * alpha;
        obs_prob(i).slice(t) = B.slice(t).cols(0, M - 1);
        obs_prob(i).slice(t).each_col() %= alpha;
      } else {
        // previous observation is missing, need to marginalize over it
        obs_prob(i).slice(t).zeros();
        for (arma::uword m = 0; m < M; m++) {
          for (arma::uword s = 0; s < S; s ++) {
            A_t.col(s) = softmax(gamma_A.slice(s) * W_A(m).slice(i).col(t));
            B_t.col(s) = softmax(gamma_B.slice(s) * W_B(m).slice(i).col(t));
          }
          Btmp = B_t.t();
          Btmp.each_col() %= A_t * alpha * arma::accu(obs_prob(i).slice(t - 1).col(m));
          obs_prob(i).slice(t) += Btmp;
        }
      }
      if(obs(t, i) < M) {
        alpha = obs_prob(i).slice(t).col(obs(t, i));
      } else {
        alpha = arma::sum(obs_prob(i).slice(t), 1);
      }
      // P(alpha_t | y_t, y_t-2, ..., y_1)
      alpha /= arma::accu(alpha);
    }
  }
}
