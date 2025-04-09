// log-likelihood and gradients of MNHMM
#include "mnhmm.h"
#include "forward.h"
#include "backward.h"
#include "eta_to_gamma.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "gradients_mnhmm.h"

// [[Rcpp::export]]
Rcpp::List log_objective_mnhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const arma::mat& X_omega,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool icpt_only_omega,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::field<arma::mat>& eta_pi,
    const arma::field<arma::cube>& eta_A,
    const Rcpp::List& eta_B,
    const arma::mat& eta_omega) {
 
  mnhmm model(
      obs, Ti, M, X_pi, X_A, X_B, X_omega, 
      icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega,
      iv_A, iv_B, tv_A, tv_B, eta_pi, eta_A, eta_B, eta_omega
  );
  arma::vec loglik(model.N);
  arma::vec loglik_i(model.D);
  arma::cube log_alpha(model.S, model.T, model.D);
  arma::cube log_beta(model.S, model.T, model.D);
  arma::mat grad_omega(model.D, model.X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(model.D);
  arma::field<arma::cube> grad_A(model.D);
  arma::field<arma::cube> grad_B(model.C, model.D);
  for (arma::uword d = 0; d < model.D; ++d) {
    grad_pi(d) = arma::mat(model.S, model.X_pi.n_rows, arma::fill::zeros);
    grad_A(d) = arma::cube(model.S, model.X_A.n_rows, model.S, arma::fill::zeros);
    for (arma::uword c = 0; c < model.C; ++c) {
      grad_B(c, d) = arma::cube(model.M(c), model.X_B(c).n_rows, model.S, arma::fill::zeros);
    }
  }
  arma::mat grad_omega2(model.D - 1, model.X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi2(model.D);
  arma::field<arma::cube> grad_A2(model.D);
  arma::field<arma::cube> grad_B2(model.C, model.D);
  for (arma::uword d = 0; d < model.D; ++d) {
    grad_pi2(d) = arma::mat(model.S - 1, model.X_pi.n_rows, arma::fill::zeros);
    grad_A2(d) = arma::cube(model.S - 1, model.X_A.n_rows, model.S, arma::fill::zeros);
    for (arma::uword c = 0; c < model.C; ++c) {
      grad_B2(c, d) = arma::cube(model.M(c) - 1, model.X_B(c).n_rows, model.S, arma::fill::zeros);
    }
  }
  arma::mat tmpmat(model.S, model.S);
  arma::mat tmpmatD(model.D, model.D);
  arma::field<arma::vec> tmpvec(model.C);
  for (arma::uword c = 0; c < model.C; ++c) {
    tmpvec(c) = arma::vec(model.M(c));
  }
  for (arma::uword i = 0; i < model.N; ++i) {
    if (!model.icpt_only_omega || i == 0) {
      model.update_omega(i);
    }
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (arma::any(model.iv_B) || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    for (arma::uword d = 0; d < model.D; ++d) {
      univariate_forward(
        log_alpha.slice(d), model.log_pi(d), model.log_A(d), 
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      univariate_backward(
        log_beta.slice(d), model.log_A(d), model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(model.Ti(i) - 1));
    }
    loglik(i) = logSumExp(model.log_omega + loglik_i);
    if (!std::isfinite(loglik(i))) {
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -model.maxval,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2),
        Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega2)
      );
    }
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (arma::uword d = 0; d < model.D; ++d) {
      // gradient wrt gamma_pi
      gradient_wrt_pi(
        grad_pi(d), tmpmat, model.log_omega, model.log_py, log_beta, loglik, 
        model.pi, model.X_pi, i, d
      );
      // gradient wrt gamma_A
      for (arma::uword t = 1; t < model.Ti(i); ++t) {
        for (arma::uword s = 0; s < model.S; ++s) {
          gradient_wrt_A(
            grad_A(d).slice(s), tmpmat, model.log_omega, model.log_py, log_alpha,
            log_beta, loglik, model.A, model.X_A, i, t, s, d
          );
        }
      }
      // gradient wrt gamma_B
      for (arma::uword c = 0; c < model.C; ++c) {
        for (arma::uword s = 0; s < model.S; ++s) {
          if (model.obs(c, 0, i) < model.M(c)) {
            gradient_wrt_B_t0(
              grad_B(c, d).slice(s), tmpvec(c), model.log_omega, model.obs, 
              model.log_pi, log_beta, loglik, model.log_B, model.B, model.X_B, 
              model.M, i, s, c, d
            );
          }
          for (arma::uword t = 1; t < model.Ti(i); ++t) {
            if (model.obs(c, t, i) < model.M(c)) {
              gradient_wrt_B(
                grad_B(c, d).slice(s), tmpvec(c), model.log_omega, model.obs,
                log_alpha, log_beta, loglik, model.log_A, model.log_B, model.B, 
                model.X_B, model.M, i, s, t, c, d
              );
            }
          }
        }
      }
    }
    gradient_wrt_omega(
      grad_omega, tmpmatD, model.omega, loglik_i, loglik, 
      model.X_omega, i
    );
  }
  grad_omega2 = model.Qd.t() * grad_omega;
  for (arma::uword d = 0; d < model.D; ++d) {
    grad_pi2(d) = model.Qs.t() * grad_pi(d);
    for (arma::uword s = 0; s < model.S; ++s) {
      grad_A2(d).slice(s) = model.Qs.t() * grad_A(d).slice(s);
    }
  }
  for (arma::uword c = 0; c < model.C; ++c) {
    for (arma::uword d = 0; d < model.D; ++d) {
      for (arma::uword s = 0; s < model.S; ++s) {
        grad_B2(c, d).slice(s) =  model.Qm(c).t() * grad_B(c, d).slice(s);
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega2)
  );
}
