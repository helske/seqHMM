// state_obs_probs algorithm for NHMM
#include "state_obs_probs.h"
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

// [[Rcpp::export]]
arma::cube predict_nhmm_singlechannel(
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
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.predict(obs_prob);
  
  return obs_prob;
}
// [[Rcpp::export]]
arma::cube predict_fanhmm_singlechannel(
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
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.predict_fanhmm(obs_prob, W_A, W_B, autoregression);
  
  return obs_prob;
}
// [[Rcpp::export]]
arma::field<arma::cube> predict_bystate_fanhmm_singlechannel(
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
arma::cube boot_predict_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::mat& gamma_pi, 
    const arma::cube& gamma_A,
    const arma::cube& gamma_B) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.predict(obs_prob);
  return obs_prob;
}
// [[Rcpp::export]]
arma::cube boot_predict_fanhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const arma::mat& gamma_pi, 
    const arma::cube& gamma_A,
    const arma::cube& gamma_B,
    const bool autoregression) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti, icpt_only_pi, icpt_only_A, 
      icpt_only_B, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  arma::cube obs_prob(model.M, model.T, model.N, arma::fill::value(arma::datum::nan));
  model.predict_fanhmm(obs_prob, W_A, W_B, autoregression);
  return obs_prob;
}
// [[Rcpp::export]]
arma::field<arma::cube> boot_predict_bystate_fanhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const arma::uvec Ti, 
    const bool icpt_only_pi, const bool icpt_only_A, const bool icpt_only_B, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const arma::mat& gamma_pi, 
    const arma::cube& gamma_A,
    const arma::cube& gamma_B,
    const bool autoregression) {
  
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

void nhmm_sc::predict(arma::cube& obs_prob) {
  
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
    obs_prob.slice(i).col(0) = B.slice(0).cols(0, M - 1).t() * alpha;
    alpha %= B.slice(0).col(obs(0, i));
    alpha /= arma::accu(alpha);
    for (arma::uword t = 1; t < Ti(i); t++) {
      alpha = A.slice(t - 1).t() * alpha;
      obs_prob.slice(i).col(t) = B.slice(t).cols(0, M - 1).t() * alpha;
      alpha %= B.slice(t).col(obs(t, i));
      alpha /= arma::accu(alpha);
    }
  }
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
      alpha = A.slice(t - 1).t() * alpha;
      obs_prob(i).slice(t) = B.slice(t).cols(0, M - 1);
      obs_prob(i).slice(t).each_col() %= alpha;
      alpha = obs_prob(i).slice(t).col(obs(t, i));
      alpha /= arma::accu(alpha);
    }
  }
}

void nhmm_sc::predict_fanhmm(
    arma::cube& obs_prob, 
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B, 
    const bool autoregression) {
  
  arma::cube A_tm1(S, S, M);
  arma::cube B_t(S, M + 1, M);
  B_t.col(M).ones(); // missing y_t
  arma::mat alpha_new(S, M); // P(alpha_t = s, y_t = m)
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
      obs_prob.slice(i).col(0).zeros();
      obs_prob(obs(0, i), 0, i) = 1.0;
    } else {
      // no autoregressive component
      obs_prob.slice(i).col(0) = B.slice(0).cols(0, M - 1).t() * alpha;
      // P(z_1) P(y_1| y_1) = P(z_1, y_1)
      alpha %= B.slice(0).col(obs(0, i));
      // P(z_1 | y_1)
      alpha /= arma::accu(alpha);
    }
    for (arma::uword t = 1; t < Ti(i); t++) {
      if (obs(t - 1, i) < M) {
        // P(alpha_t | y_t-1, ..., y_1)
        alpha = A.slice(t - 1).t() * alpha;
        // P(y_t | y_t-1,...,y_1)
        obs_prob.slice(i).col(t) = B.slice(t).cols(0, M - 1).t() * alpha;
        // P(alpha_t, y_t | y_t-1, ..., y_1) (or P(alpha_t | y_t-1, ..., y_1) if y_t missing)
        alpha %= B.slice(t).col(obs(t, i));
        // P(alpha_t | y_t, ..., y_1) (or P(alpha_t | y_t-1, ..., y_1) if y_t missing)
        alpha /= arma::accu(alpha);
      } else {
        // previous observation is missing, need to marginalize over it
        for (arma::uword m = 0; m < M; m++) {
          for (arma::uword s = 0; s < S; s ++) {
            A_tm1.slice(m).row(s) = softmax(gamma_A.slice(s) * W_A(m).slice(i).col(t - 1)).t();
            B_t.slice(m).row(s).cols(0, M - 1) = softmax(gamma_B.slice(s) * W_B(m).slice(i).col(t)).t();
          }
        }
        obs_prob.slice(i).col(t).zeros();
        for (arma::uword m = 0; m < M; m++) {
          // P(alpha_t | y_t-1 = m, y_t-2, ..., y_1)
          alpha_new.col(m) = A_tm1.slice(m).t() * alpha * obs_prob(m, t - 1, i);
          obs_prob.slice(i).col(t) += B_t.slice(m).cols(0, M - 1).t() * alpha_new.col(m);
          // P(alpha_t, y_t | y_t-1 = m, y_t-2, ..., y_1)
          alpha_new.col(m) %= B_t.slice(m).col(obs(t, i));
        }
        // P(alpha_t, y_t | y_t-2,..., y_1)
        alpha = arma::sum(alpha_new, 1);
        // P(alpha_t | y_t, y_t-2, ..., y_1)
        alpha /= arma::accu(alpha);
      }
    }
  }
}

void nhmm_sc::predict_fanhmm(
    arma::field<arma::cube>& obs_prob, 
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B, 
    const bool autoregression) {
  
  arma::cube A_tm1(S, S, M);
  arma::cube B_t(S, M + 1, M);
  B_t.col(M).ones(); // missing y_t
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
      obs_prob(i).slice(0).col(obs(0, i)) = alpha; //pi
    } else {
      // no autoregressive component
      obs_prob(i).slice(0) = B.slice(0).cols(0, M - 1);
      obs_prob(i).slice(0).each_col() %= alpha;
      alpha = obs_prob(i).slice(0).col(obs(0, i));
      alpha /= arma::accu(alpha);
    }
    for (arma::uword t = 1; t < Ti(i); t++) {
      if (obs(t - 1, i) < M) {
        // P(alpha_t | y_t-1, ..., y_1)
        alpha = A.slice(t - 1).t() * alpha;
        obs_prob(i).slice(t) = B.slice(t).cols(0, M - 1);
        obs_prob(i).slice(t).each_col() %= alpha;
        alpha = obs_prob(i).slice(t).col(obs(t, i));
        
      } else {
        // previous observation is missing, need to marginalize over it
        for (arma::uword m = 0; m < M; m++) {
          for (arma::uword s = 0; s < S; s ++) {
            A_tm1.slice(m).row(s) = softmax(gamma_A.slice(s) * W_A(m).slice(i).col(t - 1)).t();
            B_t.slice(m).row(s).cols(0, M - 1) = softmax(gamma_B.slice(s) * W_B(m).slice(i).col(t)).t();
          }
        }
        obs_prob(i).slice(t).zeros();
        for (arma::uword m = 0; m < M; m++) {
          // P(alpha_t | y_t-1 = m, y_t-2, ..., y_1)
          alpha_new.col(m) = A_tm1.slice(m).t() * alpha * arma::accu(obs_prob(i).slice(t - 1).col(m));
          Btmp = B_t.slice(m).cols(0, M - 1);
          Btmp.each_col() %= alpha_new.col(m);
          obs_prob(i).slice(t) += Btmp;
          // P(alpha_t, y_t | y_t-1 = m, y_t-2, ..., y_1)
          alpha_new.col(m) %= B_t.slice(m).col(obs(t, i));
        }
        // P(alpha_t, y_t | y_t-2,..., y_1)
        alpha = arma::sum(alpha_new, 1);
      }
      // P(alpha_t | y_t, y_t-2, ..., y_1)
      alpha /= arma::accu(alpha);
    }
  }
}
