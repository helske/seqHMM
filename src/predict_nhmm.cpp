// forward predictions for NHMMs
#include "nhmm.h"
#include "fanhmm.h"

// [[Rcpp::export]]
arma::field<arma::cube> predict_nhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B) {
  
  nhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B
  );
  arma::field<arma::cube> obs_prob(model.N, model.C);
  for (arma::uword i = 0; i < model.N; ++i) {
    for (arma::uword c = 0; c < model.C; ++c) {
      obs_prob(i, c) = arma::cube(model.S, model.M(c), model.T, arma::fill::value(arma::datum::nan));
    }
  }
  model.predict(obs_prob);
  return obs_prob;
}

// [[Rcpp::export]]
arma::field<arma::cube> boot_predict_nhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B,
    const arma::mat& gamma_pi, const arma::cube& gamma_A,
    const  arma::field<arma::cube>& gamma_B) {
  
  nhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B
  );
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  arma::field<arma::cube> obs_prob(model.N, model.C);
  for (arma::uword i = 0; i < model.N; ++i) {
    for (arma::uword c = 0; c < model.C; ++c) {
      obs_prob(i, c) = arma::cube(model.S, model.M(c), model.T, arma::fill::value(arma::datum::nan));
    }
  }
  model.predict(obs_prob);
  return obs_prob;
}

// [[Rcpp::export]]
arma::field<arma::cube> predict_fanhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const arma::field<arma::cube>& W_B1, const arma::field<arma::vec>& delta) {
  
  fanhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B,
      W_B1, delta
  );
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; ++i) {
    for (arma::uword c = 0; c < model.C; ++c) {
      obs_prob(i, c) = arma::cube(model.S, model.M(c), model.T, arma::fill::value(arma::datum::nan));
    }
  }
  model.predict(obs_prob, W_A, W_B);
  
  return obs_prob;
}
// [[Rcpp::export]]
arma::field<arma::cube> boot_predict_fanhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B,
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B,
    const arma::mat& gamma_pi, const arma::cube& gamma_A,
    const arma::field<arma::cube>& gamma_B, const arma::field<arma::cube>& W_B1, const arma::field<arma::vec>& delta) {
  
  fanhmm model(
      obs, Ti, M, X_pi, X_A, X_B, icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B,
      W_B1, delta
  );
  model.gamma_pi = gamma_pi;
  model.gamma_A = gamma_A;
  model.gamma_B = gamma_B;
  
  arma::field<arma::cube> obs_prob(model.N);
  for (arma::uword i = 0; i < model.N; ++i) {
    for (arma::uword c = 0; c < model.C; ++c) {
      obs_prob(i, c) = arma::cube(model.S, model.M(c), model.T, arma::fill::value(arma::datum::nan));
    }
  }
  model.predict(obs_prob, W_A, W_B);
  return obs_prob;
}

void nhmm::predict(arma::field<arma::cube>& obs_prob) {
  
  arma::vec alpha(S);
  arma::vec alpha_new(S);
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    alpha = pi;
    obs_prob(i, 0).slice(0) = B(0).slice(0).cols(0, M(0) - 1);
    obs_prob(i, 0).slice(0).each_col() %= alpha;
    alpha_new = obs_prob(i).slice(0).col(obs(0, 0, i));
    for (arma::uword c = 1; c < C; ++c) {
      obs_prob(i, c).slice(0) = B(c).slice(0).cols(0, M(c) - 1);
      obs_prob(i, c).slice(0).each_col() %= alpha;
      alpha_new += obs_prob(i, c).slice(0).col(obs(c, 0, i));
    }
    alpha /= arma::accu(alpha_new);
    for (arma::uword t = 1; t < Ti(i); ++t) {
      alpha = A.slice(t).t() * alpha;
      obs_prob(i, 0).slice(t) = B(0).slice(t).cols(0, M(0) - 1);
      obs_prob(i, 0).slice(t).each_col() %= alpha;
      alpha_new = obs_prob(i).slice(0).col(obs(0, t, i));
      for (arma::uword c = 1; c < C; ++c) {
        obs_prob(i, c).slice(t) = B(c).slice(t).cols(0, M(c) - 1);
        obs_prob(i, c).slice(t).each_col() %= alpha;
        alpha_new += obs_prob(i, c).slice(t).col(obs(c, t, i));
      }
      alpha /= arma::accu(alpha_new);
    }
  }
}

void fanhmm::predict(
    arma::field<arma::cube>& obs_prob, 
    const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B) {
  
  // arma::mat A_t(S, S);
  // arma::field<arma::mat> B_t(C);
  // arma::mat Btmp(S, M);
  // for (arma::uword c = 0; c < C; ++c) {
  //   B_t(c) = arma::mat(M(c), S);
  // }
  // arma::mat Btmp(S, M);
  // arma::vec alpha(S);
  // for (arma::uword i = 0; i < N; ++i) {
  //   if (!icpt_only_pi || i == 0) {
  //     update_pi(i);
  //   }
  //   if (iv_A || i == 0) {
  //     update_A(i);
  //   }
  //   if (arma::any(iv_B) || i == 0) {
  //     update_B(i);
  //   }
  //   // P(z_1)
  //   alpha = pi;
  //   // P(y_1)
  //   obs_prob(i).slice(0) = B.slice(0).cols(0, M - 1);
  //   obs_prob(i).slice(0).each_col() %= alpha;
  //   if (obs(0, i) < M) {
  //     alpha = obs_prob(i).slice(0).col(obs(0, i));
  //     alpha /= arma::accu(alpha);
  //   }
  //   for (arma::uword t = 1; t < Ti(i); ++t) {
  //     if (obs(t - 1, i) < M) {
  //       alpha = A.slice(t).t() * alpha;
  //       obs_prob(i).slice(t) = B.slice(t).cols(0, M - 1);
  //       obs_prob(i).slice(t).each_col() %= alpha;
  //     } else {
  //       // previous observation is missing, need to marginalize over it
  //       obs_prob(i).slice(t).zeros();
  //       for (arma::uword m = 0; m < M; ++m) {
  //         for (arma::uword s = 0; s < S; ++s) {
  //           A_t.col(s) = softmax(gamma_A.slice(s) * W_A(m).slice(i).col(t));
  //           B_t.col(s) = softmax(gamma_B.slice(s) * W_B(m).slice(i).col(t));
  //         }
  //         Btmp = B_t.t();
  //         Btmp.each_col() %= A_t * alpha * arma::accu(obs_prob(i).slice(t - 1).col(m));
  //         obs_prob(i).slice(t) += Btmp;
  //       }
  //     }
  //     if(obs(t, i) < M) {
  //       alpha = obs_prob(i).slice(t).col(obs(t, i));
  //     } else {
  //       alpha = arma::sum(obs_prob(i).slice(t), 1);
  //     }
  //     // P(alpha_t | y_t, y_t-2, ..., y_1)
  //     alpha /= arma::accu(alpha);
  //   }
  // }
}
