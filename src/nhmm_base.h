#ifndef NHMMBASE_H
#define NHMMBASE_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"

struct nhmm_base {
  
  const arma::uword S;
  const arma::mat& X_pi;
  const arma::cube& X_A;
  const arma::cube& X_B;
  const arma::uword K_pi;
  const arma::uword K_A;
  const arma::uword K_B;
  const arma::uword N;
  const arma::uword T;
  const arma::uvec& Ti;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const bool icpt_only_B;
  const bool iv_A;
  const bool iv_B;
  const bool tv_A; 
  const bool tv_B;
  arma::mat Qs;
  arma::mat eta_pi;
  arma::mat gamma_pi;
  arma::cube eta_A;
  arma::cube gamma_A;
  // these store pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::vec pi;
  arma::vec log_pi;
  arma::cube A;
  arma::cube log_A;
  arma::mat log_py;
  // excepted counts for EM algorithm
  arma::mat E_Pi;
  arma::field<arma::cube> E_A;
  arma::uword current_s;
  const arma::uword n_obs;
  double lambda;
  double maxval;
  int mstep_iter = 0;
  int mstep_return_code = 0;
  
  nhmm_base(
    const arma::uword S_,
    const arma::mat& X_pi_,
    const arma::cube& X_s_,
    const arma::cube& X_o_,
    const arma::uvec& Ti_, 
    const bool icpt_only_pi_,
    const bool icpt_only_A_,
    const bool icpt_only_B_,
    const bool iv_A_,
    const bool iv_B_,
    const bool tv_A_, 
    const bool tv_B_,
    arma::mat& eta_pi_,
    arma::cube& eta_A_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0,
    double maxval_ = 1e6)
    : S(S_), 
      X_pi(X_pi_),
      X_A(X_s_),
      X_B(X_o_),
      K_pi(X_pi.n_rows),
      K_A(X_A.n_rows),
      K_B(X_B.n_rows),
      N(X_A.n_slices),
      T(X_A.n_cols),
      Ti(Ti_), 
      icpt_only_pi(icpt_only_pi_),
      icpt_only_A(icpt_only_A_),
      icpt_only_B(icpt_only_B_),
      iv_A(iv_A_),
      iv_B(iv_B_),
      tv_A(tv_A_),
      tv_B(tv_B_),
      Qs(create_Q(S)),
      eta_pi(eta_pi_),
      gamma_pi(eta_to_gamma(eta_pi, Qs)), 
      eta_A(eta_A_),
      gamma_A(eta_to_gamma(eta_A, Qs)),
      pi(S),
      log_pi(S),
      A(S, S, T),
      log_A(S, S, T),
      log_py(S, T), 
      E_Pi(S, N), 
      E_A(S), 
      current_s(0),
      n_obs(n_obs_),
      lambda(lambda_),
      maxval(maxval_) {
    for (arma::uword s = 0; s < S; s++) {
      E_A(s) = arma::cube(S, N, T);
    }
  }
  void update_gamma_pi() {
    gamma_pi = eta_to_gamma(eta_pi, Qs);
  }
  void update_gamma_A() {
    gamma_A = eta_to_gamma(eta_A, Qs);
  }
  void update_pi(arma::uword i) {
    if (icpt_only_pi) {
      pi = softmax(gamma_pi.col(0));
    } else {
      pi = softmax(gamma_pi * X_pi.col(i));
    }
    log_pi = arma::log(pi);
  }
  void update_A(arma::uword i) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword s = 0; s < S; s++) { // from states
        Atmp.col(s) = softmax(gamma_A.slice(s).col(0));
      }
      A.each_slice() = Atmp.t();
    } else {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s ++) { // from states
            Atmp.col(s) = softmax(gamma_A.slice(s) * X_A.slice(i).col(t));
          }
          A.slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; s ++) { // from states
          Atmp.col(s) = softmax(gamma_A.slice(s) * X_A.slice(i).col(0));
        }
        A.each_slice() = Atmp.t();
      }
    }
    log_A = arma::log(A);
  }
  
  void estep_pi(const arma::uword i, const arma::vec& log_alpha, 
                const arma::vec& log_beta, const double ll, 
                const double pseudocount = 0) {
    E_Pi.col(i) = arma::exp(log_alpha + log_beta - ll) + pseudocount;
  }
  
  void estep_A(const arma::uword i, const arma::mat& log_alpha, 
               const arma::mat& log_beta, const double ll, 
               const double pseudocount = 0) {
    for (arma::uword k = 0; k < S; k++) { // from
      for (arma::uword j = 0; j < S; j++) { // to
        for (arma::uword t = 0; t < (Ti(i) - 1); t++) { // time
          E_A(k)(j, i, t) = exp(log_alpha(k, t) + log_A(k, j, t) + 
            log_beta(j, t + 1) + log_py(j, t + 1) - ll) + pseudocount;
        }
      }
    }
  }
  
  void mstep_pi(const double ftol_abs, const double ftol_rel, 
                const double xtol_abs, const double xtol_rel, 
                const arma::uword maxeval, const arma::uword print_level);
  void mstep_A(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const arma::uword print_level);
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
  
};

#endif
