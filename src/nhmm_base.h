#ifndef NHMMBASE_H
#define NHMMBASE_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"

struct nhmm_opt_data_pi {
  const arma::field<arma::vec>& E_Pi;
  nhmm_opt_data_pi(const arma::field<arma::vec>& E_Pi_) : E_Pi(E_Pi_) {}
};
struct nhmm_opt_data_A {
  const arma::field<arma::cube>& E_A;
  arma::uword s;
  nhmm_opt_data_A(const arma::field<arma::cube>& E_A_) : E_A(E_A_), s(0) {}
};
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
  const bool iv_pi;
  const bool iv_A;
  const bool iv_B;
  const bool tv_A; 
  const bool tv_B;
  arma::mat Qs;
  arma::mat eta_pi;
  arma::mat gamma_pi;
  arma::cube eta_A;
  arma::cube gamma_A;
  // these store Pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::vec Pi;
  arma::vec log_Pi;
  arma::cube A;
  arma::cube log_A;
  arma::mat log_py;
  
  nhmm_base(
    const arma::uword S_,
    const arma::mat& X_pi_,
    const arma::cube& X_s_,
    const arma::cube& X_o_,
    const arma::uvec& Ti_,
    const bool iv_pi_,
    const bool iv_A_,
    const bool iv_B_,
    const bool tv_A_, 
    const bool tv_B_,
    arma::mat& eta_pi_,
    arma::cube& eta_A_)
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
      iv_pi(iv_pi_),
      iv_A(iv_A_),
      iv_B(iv_B_),
      tv_A(tv_A_),
      tv_B(tv_B_),
      Qs(create_Q(S)),
      eta_pi(eta_pi_),
      gamma_pi(eta_to_gamma(eta_pi, Qs)), 
      eta_A(eta_A_),
      gamma_A(eta_to_gamma(eta_A, Qs)),
      Pi(S),
      log_Pi(S),
      A(S, S, T),
      log_A(S, S, T),
      log_py(S, T) {}
  
  void update_pi(arma::uword i) {
    Pi = softmax(gamma_pi * X_pi.col(i));
    log_Pi = arma::log(Pi);
  }
  void update_A(arma::uword i) {
    arma::mat Atmp(S, S);
    if (tv_A) {
      for (arma::uword t = 0; t < Ti(i); t++) { // time
        for (arma::uword j = 0; j < S; j ++) { // from states
          Atmp.col(j) = softmax(gamma_A.slice(j) * X_A.slice(i).col(t));
        }
        A.slice(t) = Atmp.t();
      }
    } else {
      for (arma::uword j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma_A.slice(j) * X_A.slice(i).col(0));
      }
      A.each_slice() = Atmp.t();
    }
    log_A = arma::log(A);
  }
  
  void mstep_pi(const arma::field<arma::vec>& E_Pi,
               const double xtol_abs, const double ftol_abs, const double xtol_rel,
               const double ftol_rel, arma::uword maxeval);
  void mstep_A(const arma::field<arma::cube>& E_A,
               const double xtol_abs, const double ftol_abs, const double xtol_rel,
               const double ftol_rel, arma::uword maxeval);
  double objective_pi(const arma::vec& x, arma::vec& grad, const nhmm_opt_data_pi& opt_data);
  double objective_A(const arma::vec& x, arma::vec& grad, const nhmm_opt_data_A& opt_data);
};
#endif
