// log-likelihood and gradients of NHMM
#include "nhmm_sc.h"
#include "nhmm_mc.h"
#include "mnhmm_sc.h"
#include "mnhmm_mc.h"

#include "nhmm_forward.h"
#include "nhmm_backward.h"
#include "eta_to_gamma.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "nhmm_gradients.h"

// [[Rcpp::export]]
Rcpp::List log_objective_nhmm_singlechannel(
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::cube& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const arma::uvec& Ti) {
  
  nhmm_sc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  arma::vec loglik(model.N);
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  arma::mat grad_pi(model.S - 1, model.K_pi, arma::fill::zeros);
  arma::cube grad_A(model.S - 1, model.K_A, model.S, arma::fill::zeros);
  arma::cube grad_B(model.M - 1, model.K_B, model.S, arma::fill::zeros);
  
  arma::mat tmpmat(model.S, model.S);
  arma::vec tmpvec(model.M);
  for (arma::uword i = 0; i < model.N; i++) {
    if (model.iv_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    univariate_forward_nhmm(
      log_alpha, model.log_Pi, model.log_A, 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll = logSumExp(log_alpha.col(model.Ti(i) - 1));
    if (!std::isfinite(ll)) {
      double small = -arma::datum::inf;
      grad_pi.fill(small);
      grad_A.fill(small);
      grad_B.fill(small);
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -arma::datum::inf,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B)
      );
    }
    loglik(i) = ll;
    // gradient wrt gamma_pi
    gradient_wrt_pi(grad_pi, tmpmat, model.Qs, model.log_py, log_beta, ll, model.Pi, model.X_pi, i);
    // gradient wrt gamma_A
    for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
      for (arma::uword s = 0; s < model.S; s++) {
        gradient_wrt_A(grad_A.slice(s), tmpmat, model.Qs, model.log_py, log_alpha, log_beta, ll, model.A, model.X_A, i, t, s);
      }
    }
    // gradient wrt gamma_B
    for (arma::uword s = 0; s < model.S; s++) {
      if (model.obs(0, i) < model.M) {
        gradient_wrt_B_t0(grad_B.slice(s), tmpvec, model.Qm, model.obs, model.log_Pi, log_beta, ll, model.B, model.X_B, i, s);
      }
      for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
        if (obs(t + 1, i) < model.M) {
          gradient_wrt_B(grad_B.slice(s), tmpvec, model.Qm, model.obs, log_alpha, log_beta, ll, model.log_A, model.B, model.X_B, i, s, t);
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
    arma::mat& eta_pi, const arma::mat& X_pi,
    arma::cube& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const bool iv_pi,
    const bool iv_A, const bool iv_B, const bool tv_A, const bool tv_B,
    const arma::uvec& Ti) {
  
  nhmm_mc model(
      eta_A.n_slices, X_pi, X_A, X_B, Ti,
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
  );
  
  arma::vec loglik(model.N);
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  arma::mat grad_pi(model.S - 1, model.K_pi, arma::fill::zeros);
  arma::cube grad_A(model.S - 1, model.K_A, model.S, arma::fill::zeros);
  arma::field<arma::cube> grad_B(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    grad_B(c) = arma::cube(model.M(c) - 1, model.K_B, model.S, arma::fill::zeros);
  }
  arma::mat tmpmat(model.S, model.S);
  arma::field<arma::vec> tmpvec(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    tmpvec(c) = arma::vec(model.M(c));
  }
  for (arma::uword i = 0; i < model.N; i++) {
    if (model.iv_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    univariate_forward_nhmm(
      log_alpha, model.log_Pi, model.log_A, 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward_nhmm(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll = logSumExp(log_alpha.col(model.Ti(i) - 1));
    if (!std::isfinite(ll)) {
      double small = -arma::datum::inf; // -std::max(std::min(1e10, N * 1e5), 1e3);
      grad_pi.fill(small);
      grad_A.fill(small);
      for (arma::uword c = 0; c < model.C; c++) {
        grad_B(c).fill(small);
      }
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -arma::datum::inf,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B)
      );
    }
    loglik(i) = ll;
    // gradient wrt gamma_pi
    gradient_wrt_pi(grad_pi, tmpmat, model.Qs, model.log_py, log_beta, ll, model.Pi, model.X_pi, i);
    // gradient wrt gamma_A
    for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
      for (arma::uword s = 0; s < model.S; s++) {
        gradient_wrt_A(
          grad_A.slice(s), tmpmat, model.Qs, model.log_py, log_alpha, log_beta, ll, model.A,
          model.X_A, i, t, s
        );
      }
    }
    for (arma::uword c = 0; c < model.C; c++) {
      for (arma::uword s = 0; s < model.S; s++) {
        if (model.obs(c, 0, i) < model.M(c)) {
          gradient_wrt_B_t0(
            grad_B(c).slice(s), tmpvec(c), model.Qm(c), model.obs, model.log_Pi, log_beta, ll,
            model.log_B, model.B, model.X_B, model.M, i, s, c
          );
        }
        for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
          if (model.obs(c, t + 1, i) < model.M(c)) {
            gradient_wrt_B(
              grad_B(c).slice(s), tmpvec(c), model.Qm(c), model.obs, log_alpha, log_beta,
              ll, model.log_A, model.log_B, model.B, model.X_B, model.M, i, s, t, c
            );
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
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const bool iv_omega,
    const arma::uvec& Ti) {
  
  mnhmm_sc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, iv_omega, 
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  
  arma::vec loglik(model.N);
  arma::vec loglik_i(model.D);
  arma::cube log_alpha(model.S, model.T, model.D);
  arma::cube log_beta(model.S, model.T, model.D);
  arma::mat grad_omega(model.D - 1, model.K_omega, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(model.D);
  arma::field<arma::cube> grad_A(model.D);
  arma::field<arma::cube> grad_B(model.D);
  for (arma::uword d = 0; d < model.D; d++) {
    grad_pi(d) = arma::mat(model.S - 1, model.K_pi, arma::fill::zeros);
    grad_A(d) = arma::cube(model.S - 1, model.K_A, model.S, arma::fill::zeros);
    grad_B(d) = arma::cube(model.M - 1, model.K_B, model.S, arma::fill::zeros);
  }
  arma::mat tmpmat(model.S, model.S);
  arma::mat tmpmatD(model.D, model.D);
  arma::vec tmpvec(model.M);
  for (arma::uword i = 0; i < model.N; i++) {
    if (model.iv_omega || i == 0) {
      model.update_omega(i);
    }
    if (model.iv_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    for (arma::uword d = 0; d < model.D; d++) {
      univariate_forward_nhmm(
        log_alpha.slice(d), model.log_Pi(d), model.log_A(d), 
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta.slice(d), model.log_A(d), model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(model.Ti(i) - 1));
    }
    loglik(i) = logSumExp(model.log_omega + loglik_i);
    if (!std::isfinite(loglik(i))) {
      double small = -arma::datum::inf; // -std::max(std::min(1e10, N * 1e5), 1e3);
      grad_omega.fill(small);
      for (arma::uword d = 0; d < model.D; d++) {
        grad_pi(d).fill(small);
        grad_A(d).fill(small);
        grad_B(d).fill(small);
      }
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -arma::datum::inf,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
        Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
      );
    }
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (arma::uword d = 0; d < model.D; d++) {
      // gradient wrt gamma_pi
      gradient_wrt_pi(
        grad_pi(d), tmpmat, model.Qs, model.log_omega, model.log_py, log_beta, loglik, model.Pi, model.X_pi,
        i, d
      );
      // gradient wrt gamma_A
      for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
        for (arma::uword s = 0; s < model.S; s++) {
          gradient_wrt_A(
            grad_A(d).slice(s), tmpmat, model.Qs, model.log_omega, model.log_py, log_alpha,
            log_beta, loglik, model.A, model.X_A, i, t, s, d
          );
        }
      }
      for (arma::uword s = 0; s < model.S; s++) {
        if (model.obs(0, i) < model.M) {
          gradient_wrt_B_t0(
            grad_B(d).slice(s), tmpvec, model.Qm, model.log_omega, model.obs, model.log_Pi, log_beta,
            loglik, model.B, model.X_B, i, s, d
          );
        }
        for (arma::uword t = 0; t < (model.T - 1); t++) {
          if (model.obs(t + 1, i) < model.M) {
            gradient_wrt_B(
              grad_B(d).slice(s), tmpvec, model.Qm, model.log_omega, model.obs, log_alpha,
              log_beta, loglik, model.log_A, model.B, model.X_B, i, s, t, d
            );
          }
        }
      }
    }
    gradient_wrt_omega(grad_omega, tmpmatD, model.Qd, model.omega, loglik_i, loglik, model.X_omega, i);
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
    arma::mat& eta_omega, const arma::mat& X_omega,
    arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::ucube& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
    const bool tv_A, const bool tv_B, const bool iv_omega,
    const arma::uvec& Ti) {
  
  mnhmm_mc model(
      eta_A(0).n_slices, eta_A.n_rows, X_omega, X_pi, X_A, X_B, Ti, iv_omega, 
      iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_omega, eta_pi, eta_A, eta_B
  );
  
  arma::vec loglik(model.N);
  arma::vec loglik_i(model.D);
  arma::cube log_alpha(model.S, model.T, model.D);
  arma::cube log_beta(model.S, model.T, model.D);
  arma::mat grad_omega(model.D - 1, model.K_omega, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(model.D);
  arma::field<arma::cube> grad_A(model.D);
  arma::field<arma::cube> grad_B(model.C, model.D);
  for (arma::uword d = 0; d < model.D; d++) {
    grad_pi(d) = arma::mat(model.S - 1, model.K_pi, arma::fill::zeros);
    grad_A(d) = arma::cube(model.S - 1, model.K_A, model.S, arma::fill::zeros);
    for (arma::uword c = 0; c < model.C; c++) {
      grad_B(c, d) = arma::cube(model.M(c) - 1, model.K_B, model.S, arma::fill::zeros);
    }
  }
  arma::mat tmpmat(model.S, model.S);
  arma::mat tmpmatD(model.D, model.D);
  arma::field<arma::vec> tmpvec(model.C);
  for (arma::uword c = 0; c < model.C; c++) {
    tmpvec(c) = arma::vec(model.M(c));
  }
  for (arma::uword i = 0; i < model.N; i++) {
    if (model.iv_omega || i == 0) {
      model.update_omega(i);
    }
    if (model.iv_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    for (arma::uword d = 0; d < model.D; d++) {
      univariate_forward_nhmm(
        log_alpha.slice(d), model.log_Pi(d), model.log_A(d), 
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      univariate_backward_nhmm(
        log_beta.slice(d), model.log_A(d), model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(model.Ti(i) - 1));
    }
    loglik(i) = logSumExp(model.log_omega + loglik_i);
    if (!std::isfinite(loglik(i))) {
      double small = -arma::datum::inf; // -std::max(std::min(1e10, N * 1e5), 1e3);
      grad_omega.fill(small);
      for (arma::uword d = 0; d < model.D; d++) {
        grad_pi(d).fill(small);
        grad_A(d).fill(small);
        for (arma::uword c = 0; c < model.C; c++) {
          grad_B(c, d).fill(small);
        }
      }
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -arma::datum::inf,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
        Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
      );
    }
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (arma::uword d = 0; d < model.D; d++) {
      // gradient wrt gamma_pi
      gradient_wrt_pi(
        grad_pi(d), tmpmat, model.Qs, model.log_omega, model.log_py, log_beta, loglik, model.Pi, model.X_pi,
        i, d
      );
      // gradient wrt gamma_A
      for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
        for (arma::uword s = 0; s < model.S; s++) {
          gradient_wrt_A(
            grad_A(d).slice(s), tmpmat, model.Qs, model.log_omega, model.log_py, log_alpha,
            log_beta, loglik, model.A, model.X_A, i, t, s, d
          );
        }
      }
      // gradient wrt gamma_B
      for (arma::uword c = 0; c < model.C; c++) {
        for (arma::uword s = 0; s < model.S; s++) {
          if (model.obs(c, 0, i) < model.M(c)) {
            gradient_wrt_B_t0(
              grad_B(c, d).slice(s), tmpvec(c), model.Qm(c), model.log_omega, model.obs, model.log_Pi,
              log_beta, loglik, model.log_B, model.B, model.X_B, model.M, i, s, c, d
            );
          }
          for (arma::uword t = 0; t < (model.T - 1); t++) {
            if (model.obs(c, t + 1, i) < model.M(c)) {
              gradient_wrt_B(
                grad_B(c, d).slice(s), tmpvec(c), model.Qm(c), model.log_omega, model.obs,
                log_alpha, log_beta, loglik, model.log_A, model.log_B, model.B, model.X_B, model.M, i, s, t,
                c, d
              );
            }
          }
        }
      }
    }
    gradient_wrt_omega(grad_omega, tmpmatD, model.Qd, model.omega, loglik_i, loglik, model.X_omega, i);
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega)
  );
}
