// log-likelihood and gradients of NHMM
#include "forward_nhmm.h"
#include "backward_nhmm.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "gradients.h"

// [[Rcpp::export]]
Rcpp::List log_objective_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::cube& eta_B, const arma::cube& X_o,
    const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  arma::vec loglik(N);
  arma::mat log_alpha(S, T);
  arma::mat log_beta(S, T);
  arma::mat log_py(S, T);
  
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::cube B(S, M + 1, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  
  arma::mat grad_pi(S - 1, X_i.n_rows, arma::fill::zeros);
  arma::cube grad_A(S - 1, X_s.n_rows, S, arma::fill::zeros);
  arma::cube grad_B(M - 1, X_o.n_rows, S, arma::fill::zeros);
  
  // arma::vec gradvec_S(S);
  // arma::mat gradmat_S(S, S);
  // arma::vec gradvec_M(M);
  // arma::mat gradmat_M(M, M);
  // arma::rowvec Brow(M);
  for (unsigned int i = 0; i < N; i++) {
    if (iv_pi || i == 0) {
      Pi = get_pi(eta_pi, X_i.col(i));
      log_Pi = arma::log(Pi);
    }
    if (iv_A || i == 0) {
      A = get_A(eta_A, X_s.slice(i), tv_A);
      log_A = arma::log(A);
    }
    if (iv_B || i == 0) {
      B = get_B(eta_B, X_o.slice(i), true, tv_B);
      log_B = arma::log(B);
    }
    for (unsigned int t = 0; t < Ti(i); t++) {
      log_py.col(t) = log_B.slice(t).col(obs(t, i));
    }
    log_alpha = univariate_forward_nhmm(log_Pi, log_A, log_py);
    log_beta = univariate_backward_nhmm(log_A, log_py);
    double ll = logSumExp(log_alpha.col(T - 1));
    loglik(i) = ll;
    // gradient wrt gamma_pi
    grad_pi += gradient_wrt_pi(log_py, log_beta, ll, Pi, X_i, i);
    // gradient wrt gamma_A
    for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
      for (unsigned int s = 0; s < S; s++) {
        grad_A.slice(s) += gradient_wrt_A(log_py, log_alpha, log_beta, ll, A, X_s, i, t, s);
      }
    }
    // gradient wrt gamma_B
    for (unsigned int s = 0; s < S; s++) {
      if (obs(0, i) < M) {
        grad_B.slice(s) += gradient_wrt_B_t0(obs, log_Pi, log_beta, ll, B, X_o, i, s);
      }
      for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
        if (obs(t + 1, i) < M) {
          grad_B.slice(s) += gradient_wrt_B(obs, log_alpha, log_beta, ll, log_A, B, X_o, i, s, t);
        }
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik), 
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B)
  );
}


// [[Rcpp::export]]
Rcpp::List log_objective_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::ucube& obs, const arma::uvec& M, const bool iv_pi, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uvec& Ti) {
  
  unsigned int C = M.n_elem;
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  arma::vec loglik(N);
  arma::mat log_alpha(S, T);
  arma::mat log_beta(S, T);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::mat grad_pi(S - 1, X_i.n_rows, arma::fill::zeros);
  arma::cube grad_A(S - 1, X_s.n_rows, S, arma::fill::zeros);
  arma::field<arma::cube> grad_B(C);
  for (unsigned int c = 0; c < C; c++) {
    grad_B(c) = arma::cube(M(c) - 1, X_o.n_rows, S, arma::fill::zeros);
  }
  // arma::vec gradvec_S(S);
  // arma::mat gradmat_S(S, S);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    if (iv_pi || i == 0) {
      Pi = get_pi(eta_pi, X_i.col(i));
      log_Pi = arma::log(Pi);
    }
    if (iv_A || i == 0) {
      A = get_A(eta_A, X_s.slice(i), tv_A);
      log_A = arma::log(A);
    }
    if (iv_B || i == 0) {
      B = get_B(eta_B, X_o.slice(i), M, true, tv_B);
      for (unsigned int c = 0; c < C; c++) {
        log_B(c) = arma::log(B(c));
      }
    }
    for (unsigned int t = 0; t < Ti(i); t++) {
      for (unsigned int c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
    log_alpha = univariate_forward_nhmm(log_Pi, log_A, log_py);
    log_beta = univariate_backward_nhmm(log_A, log_py);
    double ll = logSumExp(log_alpha.col(T - 1));
    loglik(i) = ll;
    // gradient wrt gamma_pi
    grad_pi += gradient_wrt_pi(log_py, log_beta, ll, Pi, X_i, i);
    // gradient wrt gamma_A
    for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
      for (unsigned int s = 0; s < S; s++) {
        grad_A.slice(s) += gradient_wrt_A(log_py, log_alpha, log_beta, ll, A, X_s, i, t, s);
      }
    }
    for (unsigned int c = 0; c < C; c++) {
      // arma::mat gradmat_M(M(c), M(c));
      // arma::rowvec Brow(M(c));
      // double logpy;
      for (unsigned int s = 0; s < S; s++) {
        if (obs(c, 0, i) < M(c)) {
          grad_B(c).slice(s) += gradient_wrt_B_t0(
            obs, log_Pi, log_beta, ll,
            log_B, B, X_o, M, i, s, c
          );
          // Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
          // gradmat_M = -Brow.t() * Brow;
          // gradmat_M.diag() += Brow;
          // logpy = 0;
          // for (unsigned int cc = 0; cc < C; cc++) {
          //   if (cc != c) {
          //     logpy += log_B(cc)(s, obs(cc, 0, i), 0);
          //   }
          // }
          // double grad = exp(log_Pi(s) + logpy + log_beta(s, 0) - ll);
          // grad_B(c).slice(s) +=
          //   gradmat_M.rows(1, M(c) - 1).col(obs(c, 0, i)) *
          //   grad * X_o.slice(i).col(0).t();
        }
        for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
          if (obs(c, t + 1, i) < M(c)) {
            grad_B(c).slice(s) += gradient_wrt_B(
              obs, log_alpha, log_beta, ll, log_A,
              log_B, B, X_o, M, i, s, t, c
            );
            // Brow = B(c).slice(t + 1).row(s).cols(0, M(c) - 1);
            // gradmat_M = -Brow.t() * Brow;
            // gradmat_M.diag() += Brow;
            // logpy = 0;
            // for (unsigned int cc = 0; cc < C; cc++) {
            //   if (cc != c) {
            //     logpy += log_B(cc)(s, obs(cc, t + 1, i), t + 1);
            //   }
            // }
            // double grad = arma::accu(
            //   exp(log_alpha.col(t) + log_A.slice(t).col(s) + 
            //     logpy + log_beta(s, t + 1) - ll));
            // grad_B(c).slice(s) += 
            //   gradmat_M.rows(1, M(c) - 1).col(obs(c, t + 1, i)) * 
            //   grad * X_o.slice(i).col(t + 1).t();
          }
        }
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B)
  );
}

// [[Rcpp::export]]
Rcpp::List log_objective_mnhmm_singlechannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const bool iv_omega, 
    const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int M = eta_B(0).n_rows + 1;
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, T, D);
  arma::cube log_beta(S, T, D);
  arma::cube log_py(S, T, D);
  
  
  arma::vec omega(D);
  arma::field<arma::vec> Pi(D);
  arma::field<arma::cube> A(D);
  arma::field<arma::cube> B(D);
  arma::vec log_omega(D);
  arma::field<arma::vec> log_Pi(D);
  arma::field<arma::cube> log_A(D);
  arma::field<arma::cube> log_B(D);
  
  arma::mat grad_omega(D - 1, X_d.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(D);
  arma::field<arma::cube> grad_A(D);
  arma::field<arma::cube> grad_B(D);
  for (unsigned int d = 0; d < D; d++) {
    grad_pi(d) = arma::mat(S - 1, X_i.n_rows, arma::fill::zeros);
    grad_A(d) = arma::cube(S - 1, X_s.n_rows, S, arma::fill::zeros);
    grad_B(d) = arma::cube(M - 1, X_o.n_rows, S, arma::fill::zeros);
  }
  arma::vec gradvec_D(D);
  arma::mat gradmat_D(D, D);
  // arma::vec gradvec_S(S);
  // arma::mat gradmat_S(S, S);
  // arma::mat gradmat_M(M + 1, M + 1);
  // arma::rowvec Brow(M);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    if (iv_omega || i == 0) {
      omega = get_omega(eta_omega, X_d.col(i));
      log_omega = arma::log(omega);
    }
    for (unsigned int d = 0; d < D; d++) {
      if (iv_pi || i == 0) {
        Pi(d) = get_pi(eta_pi(d), X_i.col(i));
        log_Pi(d) = arma::log(Pi(d));
      }
      if (iv_A || i == 0) {
        A(d) = get_A(eta_A(d), X_s.slice(i), tv_A);
        log_A(d) = arma::log(A(d));
      }
      if (iv_B || i == 0) {
        B(d) = get_B(eta_B(d), X_o.slice(i), true, tv_B);
        log_B(d) = arma::log(B(d));
      }
      for (unsigned int t = 0; t < Ti(i); t++) {
        log_py.slice(d).col(t) = log_B(d).slice(t).col(obs(t, i));
      }
      log_alpha.slice(d) = univariate_forward_nhmm(
        log_Pi(d), log_A(d), log_py.slice(d));
      log_beta.slice(d) = univariate_backward_nhmm(log_A(d), log_py.slice(d));
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(T - 1));
    }
    loglik(i) = logSumExp(log_omega + loglik_i);
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (unsigned int d = 0; d < D; d++) {
      // gradient wrt gamma_pi
      grad_pi(d) += gradient_wrt_pi(log_omega, log_py, log_beta, loglik, Pi, X_i, i, d);
      // gradient wrt gamma_A
      for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
        for (unsigned int s = 0; s < S; s++) {
          grad_A(d).slice(s) += gradient_wrt_A(log_omega, log_py, log_alpha, log_beta, loglik, A, X_s, i, t, s, d);
        }
      }
      // 
      // gradvec_S = exp(log_omega(d) + log_py.slice(d).col(0) + 
      //   log_beta.slice(d).col(0) - loglik(i));
      // // d pi / d gamma_pi
      // gradmat_S = -Pi(d) * Pi(d).t();
      // gradmat_S.diag() += Pi(d);
      // grad_pi(d) += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
      // 
      // // gradient wrt gamma_A
      // for (unsigned int t = 0; t < (T - 1); t++) {
      //   for (unsigned int s = 0; s < S; s++) {
      //     // d loglik / d a_s
      //     gradvec_S = exp(log_omega(d) + log_alpha(s, t, d) + 
      //       log_py.slice(d).col(t + 1) +
      //       log_beta.slice(d).col(t + 1) - loglik(i));
      //     // d a_s / d gamma_A
      //     gradmat_S = -A(d).slice(t).row(s).t() * A(d).slice(t).row(s);
      //     gradmat_S.diag() += A(d).slice(t).row(s);
      //     grad_A(d).slice(s) += gradmat_S.rows(1, S - 1) * 
      //       gradvec_S * X_s.slice(i).col(t).t();
      //   }
      // }
      for (unsigned int s = 0; s < S; s++) {
        if (obs(0, i) < M) {
          grad_B(d).slice(s) += gradient_wrt_B_t0(log_omega,
            obs, log_Pi, log_beta, loglik, B, X_o, i, s, d
          );
        }
        for (unsigned int t = 0; t < (T - 1); t++) {
          if (obs(t + 1, i) < M) {
            grad_B(d).slice(s) += gradient_wrt_B(log_omega,
              obs, log_alpha, log_beta, loglik, log_A, B, X_o, i, s, t, d
            );
            // gradmat_M = -B(d).slice(t + 1).row(s).cols(0, M - 1).t() * B(d).slice(t + 1).row(s).cols(0, M - 1);
            // gradmat_M.diag() += B(d).slice(t + 1).row(s).cols(0, M - 1);
            // double grad = arma::accu(
            //   exp(log_omega(d) + log_alpha.slice(d).col(t) + 
            //     log_A(d).slice(t).col(s) + log_beta(s, t + 1, d) - loglik(i)));
            // grad_B(d).slice(s) += gradmat_M.rows(1, M - 1).col(obs(t + 1, i)) * 
            //   grad * X_o.slice(i).col(t + 1).t();
          }
        }
      }
    }
    gradmat_D = -omega * omega.t();
    gradmat_D.diag() += omega;
    gradvec_D = gradmat_D * exp(loglik_i - loglik(i));
    grad_omega += (gradvec_D.rows(0, D - 2) - gradvec_D(D - 1)) * X_d.col(i).t();
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
  );
}

// [[Rcpp::export]]
Rcpp::List log_objective_mnhmm_multichannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::ucube& obs, const arma::uvec& M, const bool iv_pi, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const bool iv_omega, const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int C = M.n_elem;
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, T, D);
  arma::cube log_beta(S, T, D);
  arma::cube log_py(S, T, D);
  
  arma::vec omega(D);
  arma::field<arma::vec> Pi(D);
  arma::field<arma::cube> A(D);
  arma::field<arma::cube> B(D * C);
  arma::vec log_omega(D);
  arma::field<arma::vec> log_Pi(D);
  arma::field<arma::cube> log_A(D);
  arma::field<arma::cube> log_B(D * C);
  
  arma::mat grad_omega(D - 1, X_d.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(D);
  arma::field<arma::cube> grad_A(D);
  arma::field<arma::cube> grad_B(C, D);
  for (unsigned int d = 0; d < D; d++) {
    grad_pi(d) = arma::mat(S - 1, X_i.n_rows, arma::fill::zeros);
    grad_A(d) = arma::cube(S - 1, X_s.n_rows, S, arma::fill::zeros);
    for (unsigned int c = 0; c < C; c++) {
      grad_B(c, d) = arma::cube(M(c) - 1, X_o.n_rows, S, arma::fill::zeros);
    }
  }
  arma::vec gradvec_D(D);
  arma::mat gradmat_D(D, D);
  // arma::vec gradvec_S(S);
  // arma::mat gradmat_S(S, S);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    if (iv_omega || i == 0) {
      omega = get_omega(eta_omega, X_d.col(i));
      log_omega = arma::log(omega);
    }
    for (unsigned int d = 0; d < D; d++) {
      if (iv_pi || i == 0) {
        Pi(d) = get_pi(eta_pi(d), X_i.col(i));
        log_Pi(d) = arma::log(Pi(d));
      }
      if (iv_A || i == 0) {
        A(d) = get_A(eta_A(d), X_s.slice(i), tv_A);
        log_A(d) = arma::log(A(d));
      }
      if (iv_B || i == 0) {
        B.rows(d * C, (d + 1) * C - 1) = get_B(
          eta_B.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true, tv_B
        );
        for(unsigned int c = 0; c < C; c++) {
          log_B(d * C + c) = arma::log(B(d * C + c));
        }
      }
      for (unsigned int t = 0; t < Ti(i); t++) {
        for (unsigned int c = 0; c < C; c++) {
          log_py.slice(d).col(t) += log_B(d * C + c).slice(t).col(obs(c, t, i));
        }
      }
      log_alpha.slice(d) = univariate_forward_nhmm(
        log_Pi(d), log_A(d), log_py.slice(d)
      );
      log_beta.slice(d) = univariate_backward_nhmm(log_A(d), log_py.slice(d));
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(T - 1));
    }
    loglik(i) = logSumExp(log_omega + loglik_i);
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (unsigned int d = 0; d < D; d++) {
      grad_pi(d) += gradient_wrt_pi(log_omega, log_py, log_beta, loglik, Pi, X_i, i, d);
      // gradient wrt gamma_A
      for (unsigned int t = 0; t < (Ti(i) - 1); t++) {
        for (unsigned int s = 0; s < S; s++) {
          grad_A(d).slice(s) += gradient_wrt_A(log_omega, log_py, log_alpha, log_beta, loglik, A, X_s, i, t, s, d);
        }
      }
      // gradvec_S = exp(log_omega(d) + log_py.slice(d).col(0) + 
      //   log_beta.slice(d).col(0) - loglik(i));
      // // d pi / d gamma_pi
      // gradmat_S = -Pi(d) * Pi(d).t();
      // gradmat_S.diag() += Pi(d);
      // grad_pi(d) += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
      // 
      // // gradient wrt gamma_A
      // for (unsigned int t = 0; t < (T - 1); t++) {
      //   for (unsigned int s = 0; s < S; s++) {
      //     // d loglik / d a_s
      //     gradvec_S = exp(log_omega(d) + log_alpha(s, t, d) + 
      //       log_py.slice(d).col(t + 1) +
      //       log_beta.slice(d).col(t + 1) - loglik(i));
      //     // d a_s / d gamma_A
      //     gradmat_S = -A(d).slice(t).row(s).t() * A(d).slice(t).row(s);
      //     gradmat_S.diag() += A(d).slice(t).row(s);
      //     grad_A(d).slice(s) += gradmat_S.rows(1, S - 1) * gradvec_S * 
      //       X_s.slice(i).col(t).t();
      //   }
      // }
      // gradient wrt gamma_B
      for (unsigned int c = 0; c < C; c++) {
        // arma::mat gradmat_M(M(c), M(c));
        // arma::rowvec Brow(M(c));
        // double logpy;
        for (unsigned int s = 0; s < S; s++) {
          if (obs(c, 0, i) < M(c)) {
            grad_B(c, d).slice(s) += gradient_wrt_B_t0(
              log_omega, obs, log_Pi, log_beta, loglik, log_B, B, X_o, M, i, s, 
              c, d
            );
            // Brow = B(d * C + c).slice(0).row(s).cols(0, M(c) - 1);
            // gradmat_M = -Brow.t() * Brow;
            // gradmat_M.diag() += Brow;
            // logpy = 0;
            // for (unsigned int cc = 0; cc < C; cc++) {
            //   if (cc != c) {
            //     logpy += log_B(d * C + cc)(s, obs(cc, 0, i), 0);
            //   }
            // }
            // double grad = exp(log_omega(d) + log_Pi(d)(s) + logpy + 
            //                   log_beta(s, 0, d) - loglik(i));
            // grad_B(c, d).slice(s) += 
            //   gradmat_M.rows(1, M(c) - 1).col(obs(c, 0, i)) * 
            //   grad * X_o.slice(i).col(0).t();
          }
          for (unsigned int t = 0; t < (T - 1); t++) {
            if (obs(c, t + 1, i) < M(c)) {
              grad_B(c, d).slice(s) += gradient_wrt_B(
                log_omega, obs,log_alpha, log_beta, loglik, log_A, log_B, B, 
                X_o, M, i, s, t, c, d
              );
              // Brow = B(d * C + c).slice(t + 1).row(s).cols(0, M(c) - 1);
              // gradmat_M = -Brow.t() * Brow;
              // gradmat_M.diag() += Brow;
              // logpy = 0;
              // for (unsigned int cc = 0; cc < C; cc++) {
              //   if (cc != c) {
              //     logpy += log_B(d * C + cc)(s, obs(cc, t + 1, i), t + 1);
              //   }
              // }
              // double grad = arma::accu(
              //   exp(log_omega(d) + log_alpha.slice(d).col(t) + 
              //     log_A(d).slice(t).col(s) + logpy + log_beta(s, t + 1, d) - loglik(i)));
              // grad_B(c, d).slice(s) += 
              //   gradmat_M.rows(1, M(c) - 1).col(obs(c, t + 1, i)) * 
              //   grad * X_o.slice(i).col(t + 1).t();
            }
          }
        }
      }
    }
    gradmat_D = -omega * omega.t();
    gradmat_D.diag() += omega;
    gradvec_D = gradmat_D * exp(loglik_i - loglik(i));
    grad_omega += (gradvec_D.rows(0, D - 2) - gradvec_D(D - 1)) * X_d.col(i).t();
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
  );
}

