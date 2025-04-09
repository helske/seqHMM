// log-likelihood and gradients of NHMM
#include "nhmm.h"
#include "forward.h"
#include "backward.h"
#include "eta_to_gamma.h"
#include "get_parameters.h"
#include "logsumexp.h"
#include "gradients_nhmm.h"

// [[Rcpp::export]]
Rcpp::List log_objective_nhmm(
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
  
  arma::vec loglik(model.N);
  arma::mat log_alpha(model.S, model.T);
  arma::mat log_beta(model.S, model.T);
  
  arma::mat grad_pi(model.S, model.X_pi.n_rows, arma::fill::zeros);
  arma::cube grad_A(model.S, model.X_A.n_rows, model.S, arma::fill::zeros);
  arma::field<arma::cube> grad_B(model.C);
  for (arma::uword c = 0; c < model.C; ++c) {
    grad_B(c) = arma::cube(model.M(c), model.X_B(c).n_rows, model.S, arma::fill::zeros);
  }
  arma::mat grad_pi2(model.S - 1, model.X_pi.n_rows, arma::fill::zeros);
  arma::cube grad_A2(model.S - 1, model.X_A.n_rows, model.S, arma::fill::zeros);
  arma::field<arma::cube> grad_B2(model.C);
  for (arma::uword c = 0; c < model.C; ++c) {
    grad_B2(c) = arma::cube(model.M(c) - 1, model.X_B(c).n_rows, model.S, arma::fill::zeros);
  }
  arma::mat tmpmat(model.S, model.S);
  arma::field<arma::vec> tmpvec(model.C);
  for (arma::uword c = 0; c < model.C; ++c) {
    tmpvec(c) = arma::vec(model.M(c));
  }
  for (arma::uword i = 0; i < model.N; ++i) {
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
    univariate_forward(
      log_alpha, model.log_pi, model.log_A, 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
    univariate_backward(
      log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
    );
    double ll = logSumExp(log_alpha.col(model.Ti(i) - 1));
    if (!std::isfinite(ll)) {
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -model.maxval,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2)
      );
    }
    loglik(i) = ll;
    // gradient wrt gamma_pi
    gradient_wrt_pi(grad_pi, tmpmat, model.log_py, log_beta, ll, model.pi, model.X_pi, i);
    // gradient wrt gamma_A
    for (arma::uword t = 1; t < model.Ti(i); ++t) {
      for (arma::uword s = 0; s < model.S; ++s) {
        gradient_wrt_A(
          grad_A.slice(s), tmpmat, model.log_py, log_alpha, log_beta, ll, model.A,
          model.X_A, i, t, s
        );
      }
    }
    for (arma::uword c = 0; c < model.C; ++c) {
      for (arma::uword s = 0; s < model.S; ++s) {
        if (model.obs(c, 0, i) < model.M(c)) {
          gradient_wrt_B_t0(
            grad_B(c).slice(s), tmpvec(c), model.obs, model.log_pi, log_beta, ll,
            model.log_B, model.B, model.X_B, model.M, i, s, c
          );
        }
        for (arma::uword t = 1; t < model.Ti(i); ++t) {
          if (model.obs(c, t, i) < model.M(c)) {
            gradient_wrt_B(
              grad_B(c).slice(s), tmpvec(c), model.obs, log_alpha, log_beta,
              ll, model.log_A, model.log_B, model.B, model.X_B, model.M, i, s, t, c
            );
          }
        }
      }
    }
  }
  
  grad_pi2 = model.Qs.t() * grad_pi;
  for (arma::uword s = 0; s < model.S; ++s) {
    grad_A2.slice(s) = model.Qs.t() * grad_A.slice(s);
  }
  for (arma::uword c = 0; c < model.C; ++c) {
    for (arma::uword s = 0; s < model.S; ++s) {
      grad_B2(c).slice(s) = model.Qm(c).t() * grad_B(c).slice(s);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2)
  );
}
