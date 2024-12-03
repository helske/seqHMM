#ifndef MNHMMBASE_H
#define MNHMMBASE_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"

struct mnhmm_base {
  const arma::uword S;
  const arma::uword D;
  const arma::mat& X_omega;
  const arma::mat& X_pi;
  const arma::cube& X_A;
  const arma::cube& X_B;
  const arma::uword K_omega;
  const arma::uword K_pi;
  const arma::uword K_A;
  const arma::uword K_B;
  const arma::uword N;
  const arma::uword T;
  const arma::uvec& Ti;
  const bool icpt_only_omega;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const bool icpt_only_B;
  const bool iv_A;
  const bool iv_B;
  const bool tv_A; 
  const bool tv_B;
  arma::mat Qs;
  arma::mat Qd;
  arma::mat eta_omega;
  arma::mat gamma_omega;
  arma::field<arma::mat> eta_pi;
  arma::field<arma::mat> gamma_pi;
  arma::field<arma::cube> eta_A;
  arma::field<arma::cube> gamma_A;
  // pi, A, and log_p(y) of _one_ id we are currently working with
  arma::vec omega;
  arma::vec log_omega;
  arma::field<arma::vec> pi;
  arma::field<arma::vec> log_pi;
  arma::field<arma::cube> A;
  arma::field<arma::cube> log_A;
  arma::cube log_py;
  // excepted counts for EM algorithm
  arma::mat E_omega;
  arma::field<arma::mat> E_Pi;
  arma::field<arma::cube> E_A;
  arma::uword current_s;
  arma::uword current_d;
  const arma::uword n_obs;
  double lambda;
  double maxval;
  int mstep_iter = 0;
  int mstep_return_code = 0;
  
  mnhmm_base(
    const arma::uword S_,
    const arma::uword D_,
    const arma::mat& X_d_,
    const arma::mat& X_pi_,
    const arma::cube& X_s_,
    const arma::cube& X_o_,
    const arma::uvec& Ti_,
    const bool icpt_only_omega_,
    const bool icpt_only_pi_,
    const bool icpt_only_A_,
    const bool icpt_only_B_,
    const bool iv_A_,
    const bool iv_B_,
    const bool tv_A_, 
    const bool tv_B_,
    const arma::mat& eta_omega_,
    const arma::field<arma::mat>& eta_pi_,
    const arma::field<arma::cube>& eta_A_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0,
    double maxval_ = 1e6)
    : S(S_),
      D(D_), 
      X_omega(X_d_),
      X_pi(X_pi_),
      X_A(X_s_),
      X_B(X_o_),
      K_omega(X_omega.n_rows),
      K_pi(X_pi.n_rows),
      K_A(X_A.n_rows),
      K_B(X_B.n_rows),
      N(X_A.n_slices),
      T(X_A.n_cols),
      Ti(Ti_),
      icpt_only_omega(icpt_only_omega_),
      icpt_only_pi(icpt_only_pi_),
      icpt_only_A(icpt_only_A_),
      icpt_only_B(icpt_only_B_),
      iv_A(iv_A_),
      iv_B(iv_B_),
      tv_A(tv_A_),
      tv_B(tv_B_),
      Qs(create_Q(S)),
      Qd(create_Q(D)),
      eta_omega(eta_omega_),
      gamma_omega(eta_to_gamma(eta_omega, Qd)),
      eta_pi(eta_pi_),
      gamma_pi(eta_to_gamma(eta_pi, Qs)), 
      eta_A(eta_A_),
      gamma_A(eta_to_gamma(eta_A, Qs)),
      omega(D),
      log_omega(D),
      pi(D),
      log_pi(D),
      A(D),
      log_A(D),
      log_py(S, T, D), 
      E_omega(D, N),
      E_Pi(D),
      E_A(S, D),
      current_s(0),
      current_d(0),
      n_obs(n_obs_),
      lambda(lambda_),
      maxval(maxval_) {
    for (arma::uword d = 0; d < D; d++) {
      pi(d) = arma::vec(S);
      log_pi(d) = arma::vec(S);
      A(d) = arma::cube(S, S, T);
      log_A(d) = arma::cube(S, S, T);
      E_Pi(d) = arma::mat(S, N);
      for (arma::uword s = 0; s < S; s++) {
        E_A(s, d) = arma::cube(S, N, T);
      }
    }
  }
  void update_gamma_omega() {
    gamma_omega = eta_to_gamma(eta_omega, Qd);
    
  }
  void update_gamma_pi() {
    for (arma::uword d = 0; d < D; d++) {
      gamma_pi(d) = eta_to_gamma(eta_pi(d), Qs);
    }
  }
  void update_gamma_A() {
    for (arma::uword d = 0; d < D; d++) {
      gamma_A(d) = eta_to_gamma(eta_A(d), Qs);
    }
  }
  void update_omega(const arma::uword i) {
    if (icpt_only_omega) {
      omega = softmax(gamma_omega.col(0));
    } else {
      omega = softmax(gamma_omega * X_omega.col(i));  
    }
    log_omega = arma::log(omega);
  }
  void update_pi(const arma::uword i) {
    if (icpt_only_pi) {
      for (arma::uword d = 0; d < D; d++) {
        pi(d) = softmax(gamma_pi(d).col(0));
        log_pi(d) = arma::log(pi(d));
      }
    } else {
      for (arma::uword d = 0; d < D; d++) {
        pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
        log_pi(d) = arma::log(pi(d));
      }
    }
  }
  void update_pi(const arma::uword i, const arma::uword d) {
    if (icpt_only_pi) {
      pi(d) = softmax(gamma_pi(d).col(0));
    } else {
      pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
    }
    log_pi(d) = arma::log(pi(d));
  }
  void update_A(const arma::uword i) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword d = 0; d < D; d++) {
        for (arma::uword s = 0; s < S; s++) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s).col(0));
        }
        A(d).each_slice() = Atmp.t();
        log_A(d) = arma::log(A(d));
      }
    } else {
      for (arma::uword d = 0; d < D; d++) {
        if (tv_A) {
          for (arma::uword t = 0; t < Ti(i); t++) { // time
            for (arma::uword s = 0; s < S; s++) { // from states
              Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(t));
            }
            A(d).slice(t) = Atmp.t();
          }
        } else {
          for (arma::uword s = 0; s < S; s++) { // from states
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(0));
          }
          A(d).each_slice() = Atmp.t();
        }
        log_A(d) = arma::log(A(d));
      }
    }
  }
  void update_A(const arma::uword i, const arma::uword d) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword s = 0; s < S; s++) { // from states
        Atmp.col(s) = softmax(gamma_A(d).slice(s).col(0));
      }
      A(d).each_slice() = Atmp.t();
    } else {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s++) { // from states
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(t));
          }
          A(d).slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; s++) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(0));
        }
        A(d).each_slice() = Atmp.t();
      }
      
    }
    log_A(d) = arma::log(A(d));
  }
  
  void estep_omega(const arma::uword i, const arma::vec ll_i, 
                   const double ll) {
    E_omega.col(i) = arma::exp(ll_i - ll);
  }
  
  void estep_pi(const arma::uword i, const arma::uword d, 
                const arma::vec& log_alpha, 
                const arma::vec& log_beta, const double ll) {
    E_Pi(d).col(i) = arma::exp(log_alpha + log_beta - ll);
  }
  
  void estep_A(const arma::uword i, const arma::uword d, 
               const arma::mat& log_alpha, const arma::mat& log_beta, 
               const double ll) {
    for (arma::uword k = 0; k < S; k++) { // from
      for (arma::uword j = 0; j < S; j++) { // to
        for (arma::uword t = 0; t < (Ti(i) - 1); t++) { // time
          E_A(k, d)(j, i, t) = exp(log_alpha(k, t) + log_A(d)(k, j, t) + 
            log_beta(j, t + 1) + log_py(j, t + 1, d) - ll);
        }
      }
    }
  }
  void mstep_omega(const double ftol_abs, const double ftol_rel, 
                   const double xtol_abs, const double xtol_rel, 
                   const arma::uword maxeval, const double bound, 
                   const arma::uword print_level);
  void mstep_pi(const double ftol_abs, const double ftol_rel, 
                const double xtol_abs, const double xtol_rel, 
                const arma::uword maxeval, const double bound, 
                const arma::uword print_level);
  void mstep_A(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const double bound, 
               const arma::uword print_level);
  double objective_omega(const arma::vec& x, arma::vec& grad);
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
};
#endif
