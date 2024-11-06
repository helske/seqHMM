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
  const bool iv_omega;
  const bool iv_pi;
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
  // Pi, A, and log_p(y) of _one_ id and cluster we are currently working with
  arma::vec omega;
  arma::vec log_omega;
  arma::field<arma::vec> Pi;
  arma::field<arma::vec> log_Pi;
  arma::field<arma::cube> A;
  arma::field<arma::cube> log_A;
  arma::cube log_py;
  arma::uword n_obs;
  double penalty;
  mnhmm_base(
    const arma::uword S_,
    const arma::uword D_,
    const arma::mat& X_d_,
    const arma::mat& X_pi_,
    const arma::cube& X_s_,
    const arma::cube& X_o_,
    const arma::uvec& Ti_,
    const bool iv_omega_,
    const bool iv_pi_,
    const bool iv_A_,
    const bool iv_B_,
    const bool tv_A_, 
    const bool tv_B_,
    arma::mat& eta_omega_,
    arma::field<arma::mat>& eta_pi_,
    arma::field<arma::cube>& eta_A_,
    const double penalty = 0)
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
      iv_omega(iv_omega_),
      iv_pi(iv_pi_),
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
      Pi(D),
      log_Pi(D),
      A(D),
      log_A(D),
      log_py(S, T, D), 
      n_obs(sum(Ti)),
      penalty(penalty){
    for (arma::uword d = 0; d < D; d++) {
      Pi(d) = arma::vec(S);
      log_Pi(d) = arma::vec(S);
      A(d) = arma::cube(S, S, T);
      log_A(d) = arma::cube(S, S, T);
    }
  }
  
  void update_omega(arma::uword i) {
    omega = softmax(gamma_omega * X_omega.col(i));
    log_omega = arma::log(omega);
  }
  void update_pi(arma::uword i) {
    for (arma::uword d = 0; d < D; d++) {
      Pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
      log_Pi(d) = arma::log(Pi(d));
    }
  }
  void update_A(arma::uword i) {
    arma::mat Atmp(S, S);
    for (arma::uword d = 0; d < D; d++) {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword j = 0; j < S; j ++) { // from states
            Atmp.col(j) = softmax(gamma_A(d).slice(j) * X_A.slice(i).col(t));
          }
          A(d).slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword j = 0; j < S; j ++) { // from states
          Atmp.col(j) = softmax(gamma_A(d).slice(j) * X_A.slice(i).col(0));
        }
        A(d).each_slice() = Atmp.t();
      }
      log_A(d) = arma::log(A(d));
    }
  }
};
#endif
