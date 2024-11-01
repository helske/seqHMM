#ifndef NHMM_H
#define NHMM_H

#include <RcppArmadillo.h>
#include "eta_to_gamma.h"

struct nhmm {
  const arma::umat& obs;
  const arma::mat& X_pi;
  const arma::cube& X_A;
  const arma::cube& X_B;  
  const arma::mat& Qs;
  const arma::mat& Qm;
  const arma::uvec& Ti;
  arma::mat eta_pi;
  arma::cube eta_A;
  arma::cube eta_B;
  const arma::uword N;
  const arma::uword T;
  const arma::uword S;
  const arma::uword M;
  const bool iv_pi;
  const bool iv_A;
  const bool iv_B;
  const bool tv_A; 
  const bool tv_B;
  const arma::uword np_pi;
  const arma::uword np_A;
  const arma::uword np_B;
  const arma::uword K_pi;
  const arma::uword K_A;
  const arma::uword K_B;
  arma::mat gamma_pi;
  arma::cube gamma_A;
  arma::cube gamma_B;
  
  nhmm(const arma::umat& obs_,
       const arma::mat& X_pi_,
       const arma::cube& X_s_,
       const arma::cube& X_o_,
       const arma::mat& Qs_,
       const arma::mat& Qm_,
       const arma::uvec& Ti_,
       arma::mat& eta_pi_,
       arma::cube& eta_A_,
       arma::cube& eta_B_,
       const bool iv_pi_,
       const bool iv_A_,
       const bool iv_B_,
       const bool tv_A_, 
       const bool tv_B_)
    : obs(obs_),
      X_pi(X_pi_),
      X_A(X_s_),
      X_B(X_o_),
      Qs(Qs_),
      Qm(Qm_),
      Ti(Ti_),
      eta_pi(eta_pi_),
      eta_A(eta_A_),
      eta_B(eta_B_),
      N(obs.n_cols),
      T(obs.n_rows),
      S(Qs.n_rows),
      M(Qm.n_rows),
      iv_pi(iv_pi_),
      iv_A(iv_A_),
      iv_B(iv_B_),
      tv_A(tv_A_),
      tv_B(tv_B_),
      np_pi(eta_pi.n_elem),
      np_A(eta_A.n_elem),
      np_B(eta_B.n_elem),
      K_pi(X_pi.n_rows),
      K_A(X_A.n_rows),
      K_B(X_B.n_rows) {
    gamma_pi = eta_to_gamma(eta_pi, Qs);
    gamma_A = eta_to_gamma(eta_A, Qs);
    gamma_B = eta_to_gamma(eta_B, Qs);
  }
  void mstep(
      const arma::vec E_Pi, const arma::cube E_A, const arma::cube E_B,
      const double xtol_abs, const double ftol_abs, const double xtol_rel, 
      const double ftol_rel, arma::uword maxeval);
};
#endif
