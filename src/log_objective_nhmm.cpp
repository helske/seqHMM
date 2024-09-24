// log-likelihood and gradients of NHMM
#include "forward_nhmm.h"
#include "backward_nhmm.h"
#include "get_parameters.h"
#include "logsumexp.h"

// [[Rcpp::export]]
Rcpp::List log_objective_nhmm_singlechannel(
    const arma::mat& gamma_pi_raw, const arma::mat& X_i,
    const arma::cube& gamma_A_raw, const arma::cube& X_s,
    const arma::cube& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw.n_slices;
  unsigned int M = gamma_B_raw.n_rows + 1;
  arma::vec loglik(N);
  arma::mat log_alpha(S, T);
  arma::mat log_beta(S, T);
  arma::mat log_py(S, T);
  
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  
  arma::mat grad_pi(S - 1, X_i.n_rows, arma::fill::zeros);
  arma::cube grad_A(S - 1, X_s.n_rows, S, arma::fill::zeros);
  arma::cube grad_B(M - 1, X_o.n_rows, S, arma::fill::zeros);
  
  arma::vec gradvec_S(S);
  arma::mat gradmat_S(S, S);
  arma::mat gradmat_M(M, M);
  arma::mat A(S, S);
  arma::rowvec Brow(M);
  for (unsigned int i = 0; i < N; i++) {
    if (iv_pi || i == 0) {
      log_Pi = get_pi(gamma_pi_raw, X_i.col(i), true);
    }
    if (iv_A || i == 0) {
      log_A = get_A(gamma_A_raw, X_s.slice(i), true, tv_A);
    }
    if (iv_B || i == 0) {
      log_B = get_B(gamma_B_raw, X_o.slice(i), true, true, tv_B);
    }
    for (unsigned int t = 0; t < Ti[i]; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(s, obs(t, i), t);
      }
    }
    log_alpha = univariate_forward_nhmm(log_Pi, log_A, log_py);
    log_beta = univariate_backward_nhmm(log_A, log_py);
    double ll = logSumExp(log_alpha.col(T - 1));
    loglik(i) = ll;
    // gradient wrt gamma_pi
    // d loglik / d pi
    gradvec_S = exp(log_py.col(0) + log_beta.col(0) - ll);
    // d pi / d gamma_pi
    arma::vec Pi = exp(log_Pi);
    gradmat_S = -Pi * Pi.t();
    gradmat_S.diag() += Pi;
    //gradmat_S = arma::diagmat(Pi) - Pi * Pi.t();
    grad_pi += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
    
    // gradient wrt gamma_A
    for (unsigned int t = 0; t < (Ti[i] - 1); t++) {
      A = exp(log_A.slice(t));
      for (unsigned int s = 0; s < S; s++) {
        // d loglik / d a_s
        gradvec_S = exp(log_alpha(s, t) + log_py.col(t + 1) + log_beta.col(t + 1) - ll);
        // d a_s / d gamma_A
        gradmat_S = -A.row(s).t() * A.row(s);
        gradmat_S.diag() += A.row(s);
        //arma::diagmat(A.row(s)) - A.row(s).t() * A.row(s);
        grad_A.slice(s) += gradmat_S.rows(1, S - 1) * gradvec_S * X_s.slice(i).col(t).t();
      }
    }
    // gradient wrt gamma_B
    for (unsigned int s = 0; s < S; s++) {
      if (obs(0, i) < M) {
        Brow = exp(log_B.slice(0).row(s).cols(0, M - 1));
        gradmat_M = -Brow.t() * Brow;
        gradmat_M.diag() += Brow;
        double grad = exp(log_Pi(s) + log_beta(s, 0) - ll);
        grad_B.slice(s) += gradmat_M.rows(1, M - 1).col(obs(0, i)) * grad * X_o.slice(i).col(0).t();
      }
      for (unsigned int t = 0; t < (Ti[i] - 1); t++) {
        if (obs(t + 1, i) < M) {
          Brow = exp(log_B.slice(t + 1).row(s).cols(0, M - 1));
          gradmat_M = -Brow.t() * Brow;
          gradmat_M.diag() += Brow;
          double grad = arma::accu(
            exp(log_alpha.col(t) + log_A.slice(t).col(s) + log_beta(s, t + 1) - ll));
          grad_B.slice(s) += gradmat_M.rows(1, M - 1).col(obs(t + 1, i)) * grad * X_o.slice(i).col(t + 1).t();
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
    const arma::mat& gamma_pi_raw, const arma::mat& X_i,
    const arma::cube& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::cube& obs, const arma::uvec& M, const bool iv_pi, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uvec& Ti) {
  
  unsigned int C = M.n_elem;
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw.n_slices;
  arma::vec loglik(N);
  arma::mat log_alpha(S, T);
  arma::mat log_beta(S, T);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::mat grad_pi(S - 1, X_i.n_rows, arma::fill::zeros);
  arma::cube grad_A(S - 1, X_s.n_rows, S, arma::fill::zeros);
  arma::field<arma::cube> grad_B(C);
  for (unsigned int c = 0; c < C; c++) {
    grad_B(c) = arma::cube(M(c) - 1, X_o.n_rows, S, arma::fill::zeros);
  }
  arma::vec gradvec_S(S);
  arma::mat gradmat_S(S, S);
  arma::mat A(S, S);
  for (unsigned int i = 0; i < N; i++) {
    if (iv_pi || i == 0) {
      log_Pi = get_pi(gamma_pi_raw, X_i.col(i), true);
    }
    if (iv_A || i == 0) {
      log_A = get_A(gamma_A_raw, X_s.slice(i), true, tv_A);
    }
    if (iv_B || i == 0) {
      log_B = get_B(gamma_B_raw, X_o.slice(i), M, true, true, tv_B);
    }
    for (unsigned int t = 0; t < Ti[i]; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = 0;
        for (unsigned int c = 0; c < C; c++) {
          log_py(s, t) += log_B(c)(s, obs(c, t, i), t);
        }
      }
    }
    log_alpha = univariate_forward_nhmm(log_Pi, log_A, log_py);
    log_beta = univariate_backward_nhmm(log_A, log_py);
    double ll = logSumExp(log_alpha.col(T - 1));
    loglik(i) = ll;
    // gradient wrt gamma_pi
    // d loglik / d pi
    gradvec_S = exp(log_py.col(0) + log_beta.col(0) - ll);
    // d pi / d gamma_pi
    arma::vec Pi = exp(log_Pi);
    gradmat_S = -Pi * Pi.t();
    gradmat_S.diag() += Pi;
    grad_pi += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
    // gradient wrt gamma_A
    for (unsigned int t = 0; t < (Ti[i] - 1); t++) {
      A = exp(log_A.slice(t));
      for (unsigned int s = 0; s < S; s++) {
        // d loglik / d a_s
        gradvec_S = exp(log_alpha(s, t) + log_py.col(t + 1) + 
          log_beta.col(t + 1) - ll);
        // d a_s / d gamma_A
        gradmat_S = -A.row(s).t() * A.row(s);
        gradmat_S.diag() += A.row(s);
        grad_A.slice(s) += gradmat_S.rows(1, S - 1) * 
          gradvec_S * X_s.slice(i).col(t).t();
      }
    }
    for (unsigned int c = 0; c < C; c++) {
      arma::mat gradmat_M(M(c), M(c));
      arma::rowvec Brow(M(c));
      double logpy;
      for (unsigned int s = 0; s < S; s++) {
        if (obs(c, 0, i) < M(c)) {
          Brow = exp(log_B(c).slice(0).row(s).cols(0, M(c) - 1));
          gradmat_M = -Brow.t() * Brow;
          gradmat_M.diag() += Brow;
          logpy = 0;
          for (unsigned int cc = 0; cc < C; cc++) {
            if (cc != c) {
              logpy += log_B(cc)(s, obs(cc, 0, i), 0);
            }
          }
          double grad = exp(log_Pi(s) + logpy + log_beta(s, 0) - ll);
          grad_B(c).slice(s) += 
            gradmat_M.rows(1, M(c) - 1).col(obs(c, 0, i)) * 
            grad * X_o.slice(i).col(0).t();
        }
        for (unsigned int t = 0; t < (Ti[i] - 1); t++) {
          if (obs(c, t + 1, i) < M(c)) {
            Brow = exp(log_B(c).slice(t + 1).row(s).cols(0, M(c) - 1));
            gradmat_M = -Brow.t() * Brow;
            gradmat_M.diag() += Brow;
            logpy = 0;
            for (unsigned int cc = 0; cc < C; cc++) {
              if (cc != c) {
                logpy += log_B(cc)(s, obs(cc, t + 1, i), t + 1);
              }
            }
            double grad = arma::accu(
              exp(log_alpha.col(t) + log_A.slice(t).col(s) + 
                logpy + log_beta(s, t + 1) - ll));
            grad_B(c).slice(s) += 
              gradmat_M.rows(1, M(c) - 1).col(obs(c, t + 1, i)) * 
              grad * X_o.slice(i).col(t + 1).t();
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
    const arma::field<arma::mat>& gamma_pi_raw, const arma::mat& X_i,
    const arma::field<arma::cube>& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& gamma_omega_raw, const arma::mat& X_d,
    const arma::mat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const bool iv_omega, 
    const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw(0).n_slices;
  unsigned int D = gamma_omega_raw.n_rows + 1;
  unsigned int M = gamma_B_raw(0).n_rows + 1;
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, T, D);
  arma::cube log_beta(S, T, D);
  arma::cube log_py(S, T, D);
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
  arma::vec gradvec_S(S);
  arma::mat gradmat_S(S, S);
  arma::mat gradmat_M(M + 1, M + 1);
  arma::vec Pi(S);
  arma::mat A(S, S);
  arma::rowvec Brow(M);
  arma::vec omega(D);
  for (unsigned int i = 0; i < N; i++) {
    if (iv_omega || i == 0) {
      log_omega = get_omega(gamma_omega_raw, X_d.col(i), true);
    }
    for (unsigned int d = 0; d < D; d++) {
      if (iv_pi || i == 0) {
        log_Pi(d) = get_pi(gamma_pi_raw(d), X_i.col(i), true);
      }
      if (iv_A || i == 0) {
        log_A(d) = get_A(gamma_A_raw(d), X_s.slice(i), true, tv_A);
      }
      if (iv_B || i == 0) {
        log_B(d) = get_B(gamma_B_raw(d), X_o.slice(i), true, true, tv_B);
      }
      for (unsigned int t = 0; t < Ti[i]; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(s, t, d) = log_B(d)(s, obs(t, i), t);
        }
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
      gradvec_S = exp(log_omega(d) + log_py.slice(d).col(0) + 
        log_beta.slice(d).col(0) - loglik(i));
      // d pi / d gamma_pi
      Pi = exp(log_Pi(d));
      gradmat_S = -Pi * Pi.t();
      gradmat_S.diag() += Pi;
      grad_pi(d) += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
      
      // gradient wrt gamma_A
      for (unsigned int t = 0; t < (T - 1); t++) {
        A = exp(log_A(d).slice(t));
        for (unsigned int s = 0; s < S; s++) {
          // d loglik / d a_s
          gradvec_S = exp(log_omega(d) + log_alpha(s, t, d) + 
            log_py.slice(d).col(t + 1) +
            log_beta.slice(d).col(t + 1) - loglik(i));
          // d a_s / d gamma_A
          gradmat_S = -A.row(s).t() * A.row(s);
          gradmat_S.diag() += A.row(s);
          grad_A(d).slice(s) += gradmat_S.rows(1, S - 1) * 
            gradvec_S * X_s.slice(i).col(t).t();
        }
      }
      for (unsigned int s = 0; s < S; s++) {
        if (obs(0, i) < M) {
          Brow = exp(log_B(d).slice(0).row(s).cols(0, M - 1));
          gradmat_M = -Brow.t() * Brow;
          gradmat_M.diag() += Brow;
          double grad = exp(log_omega(d) + log_Pi(d)(s) + 
                            log_beta(s, 0, d) - loglik(i));
          grad_B(d).slice(s) += gradmat_M.rows(1, M - 1).col(obs(0, i)) * 
            grad * X_o.slice(i).col(0).t();
        }
        for (unsigned int t = 0; t < (T - 1); t++) {
          if (obs(t + 1, i) < M) {
            Brow = exp(log_B(d).slice(t + 1).row(s).cols(0, M - 1));
            gradmat_M = -Brow.t() * Brow;
            gradmat_M.diag() += Brow;
            double grad = arma::accu(
              exp(log_omega(d) + log_alpha.slice(d).col(t) + 
                log_A(d).slice(t).col(s) + log_beta(s, t + 1, d) - loglik(i)));
            grad_B(d).slice(s) += gradmat_M.rows(1, M - 1).col(obs(t + 1, i)) * 
              grad * X_o.slice(i).col(t + 1).t();
          }
        }
      }
    }
    gradvec_D = exp(loglik_i - loglik(i));
    omega = exp(log_omega);
    gradmat_D = -omega * omega.t();
    gradmat_D.diag() += omega;
    grad_omega += gradmat_D.rows(1, D - 1) * gradvec_D * X_d.col(i).t();
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
    const arma::field<arma::mat>& gamma_pi_raw, const arma::mat& X_i,
    const arma::field<arma::cube>& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& gamma_omega_raw, const arma::mat& X_d,
    const arma::cube& obs, const arma::uvec& M, const bool iv_pi, 
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const bool iv_omega, const arma::uvec& Ti) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw(0).n_slices;
  unsigned int D = gamma_omega_raw.n_rows + 1;
  unsigned int C = M.n_elem;
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, T, D);
  arma::cube log_beta(S, T, D);
  arma::cube log_py(S, T, D);
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
  arma::vec gradvec_S(S);
  arma::mat gradmat_S(S, S);
  arma::vec Pi(S);
  arma::mat A(S, S);
  arma::vec omega(D);
  for (unsigned int i = 0; i < N; i++) {
    if (iv_omega || i == 0) {
      log_omega = get_omega(gamma_omega_raw, X_d.col(i), true);
    }
    for (unsigned int d = 0; d < D; d++) {
      if (iv_pi || i == 0) {
        log_Pi(d) = get_pi(gamma_pi_raw(d), X_i.col(i), true);
      }
      if (iv_A || i == 0) {
        log_A(d) = get_A(gamma_A_raw(d), X_s.slice(i), true, tv_A);
      }
      if (iv_B || i == 0) {
        log_B.rows(d * C, (d + 1) * C - 1) = get_B(
          gamma_B_raw.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true, 
          true, tv_B
        );
      }
      for (unsigned int t = 0; t < Ti[i]; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(s, t, d) = 0;
          for (unsigned int c = 0; c < C; c++) {
            log_py(s, t, d) += log_B(d * C + c)(s, obs(c, t, i), t);
          }
        }
      }
      log_alpha.slice(d) = univariate_forward_nhmm(log_Pi(d), log_A(d), 
                      log_py.slice(d));
      log_beta.slice(d) = univariate_backward_nhmm(log_A(d), log_py.slice(d));
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(T - 1));
    }
    loglik(i) = logSumExp(log_omega + loglik_i);
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (unsigned int d = 0; d < D; d++) {
      gradvec_S = exp(log_omega(d) + log_py.slice(d).col(0) + 
        log_beta.slice(d).col(0) - loglik(i));
      // d pi / d gamma_pi
      Pi = exp(log_Pi(d));
      gradmat_S = -Pi * Pi.t();
      gradmat_S.diag() += Pi;
      grad_pi(d) += gradmat_S.rows(1, S - 1) * gradvec_S * X_i.col(i).t();
      
      // gradient wrt gamma_A
      for (unsigned int t = 0; t < (T - 1); t++) {
        A = exp(log_A(d).slice(t));
        for (unsigned int s = 0; s < S; s++) {
          // d loglik / d a_s
          gradvec_S = exp(log_omega(d) + log_alpha(s, t, d) + 
            log_py.slice(d).col(t + 1) +
            log_beta.slice(d).col(t + 1) - loglik(i));
          // d a_s / d gamma_A
          gradmat_S = -A.row(s).t() * A.row(s);
          gradmat_S.diag() += A.row(s);
          grad_A(d).slice(s) += gradmat_S.rows(1, S - 1) * gradvec_S * 
            X_s.slice(i).col(t).t();
        }
      }
      // gradient wrt gamma_B
      for (unsigned int c = 0; c < C; c++) {
        arma::mat gradmat_M(M(c), M(c));
        arma::rowvec Brow(M(c));
        double logpy;
        for (unsigned int s = 0; s < S; s++) {
          if (obs(c, 0, i) < M(c)) {
            Brow = exp(log_B(d * C + c).slice(0).row(s).cols(0, M(c) - 1));
            gradmat_M = -Brow.t() * Brow;
            gradmat_M.diag() += Brow;
            logpy = 0;
            for (unsigned int cc = 0; cc < C; cc++) {
              if (cc != c) {
                logpy += log_B(d * C + cc)(s, obs(cc, 0, i), 0);
              }
            }
            double grad = exp(log_omega(d) + log_Pi(d)(s) + logpy + 
                              log_beta(s, 0, d) - loglik(i));
            grad_B(c, d).slice(s) += 
              gradmat_M.rows(1, M(c) - 1).col(obs(c, 0, i)) * 
              grad * X_o.slice(i).col(0).t();
          }
          for (unsigned int t = 0; t < (T - 1); t++) {
            if (obs(c, t + 1, i) < M(c)) {
              Brow = exp(log_B(d * C + c).slice(t + 1).row(s).cols(0, M(c) - 1));
              gradmat_M = -Brow.t() * Brow;
              gradmat_M.diag() += Brow;
              logpy = 0;
              for (unsigned int cc = 0; cc < C; cc++) {
                if (cc != c) {
                  logpy += log_B(d * C + cc)(s, obs(cc, t + 1, i), t + 1);
                }
              }
              double grad = arma::accu(
                exp(log_omega(d) + log_alpha.slice(d).col(t) + 
                  log_A(d).slice(t).col(s) + logpy + log_beta(s, t + 1, d) - loglik(i)));
              grad_B(c, d).slice(s) += 
                gradmat_M.rows(1, M(c) - 1).col(obs(c, t + 1, i)) * 
                grad * X_o.slice(i).col(t + 1).t();
            }
          }
        }
      }
    }
    gradvec_D = exp(loglik_i - loglik(i));
    omega = exp(log_omega);
    gradmat_D = -omega * omega.t();
    gradmat_D.diag() += omega;
    grad_omega += gradmat_D.rows(1, D - 1) * gradvec_D * X_d.col(i).t();
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
  );
}

