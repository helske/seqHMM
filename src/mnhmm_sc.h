#ifndef MNHMM_SC_H
#define MNHMM_SC_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "mnhmm_base.h"

struct mnhmm_sc : public mnhmm_base {
  
  const arma::umat& obs;
  arma::field<arma::cube> eta_B;
  const arma::uword M;
  arma::mat Qm;
  arma::field<arma::cube> gamma_B;
  // these store Pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  
  mnhmm_sc(
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
    const arma::umat& obs_,
    arma::mat& eta_omega_,
    arma::field<arma::mat>& eta_pi_,
    arma::field<arma::cube>& eta_A_,
    arma::field<arma::cube>& eta_B_,
    const double penalty = 0)
    : mnhmm_base(S_, D_, X_d_, X_pi_, X_s_, X_o_, Ti_, iv_omega_, iv_pi_, iv_A_, 
      iv_B_, tv_A_, tv_B_, eta_omega_, eta_pi_, eta_A_, penalty),
      obs(obs_), 
      eta_B(eta_B_), 
      M(eta_B(0).n_rows + 1),
      Qm(create_Q(M)), 
      gamma_B(eta_to_gamma(eta_B, Qm)),
      B(D),
      log_B(D){
    for (arma::uword d = 0; d < D; d++) {
      B(d) = arma::cube(S, M + 1, T);
      log_B(d) = arma::cube(S, M + 1, T);
    }
  }
  
  void update_B(const arma::uword i) {
    arma::mat Btmp(M + 1, S, arma::fill::ones);
    if (tv_B) {
      for (arma::uword d = 0; d < D; d++) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s++) { // from states
            Btmp.col(s).rows(0, M - 1) = softmax(gamma_B(d).slice(s) * X_B.slice(i).col(t));
          }
          B(d).slice(t) = Btmp.t();
        }
        log_B(d) = arma::log(B(d));
      }
    } else {
      for (arma::uword d = 0; d < D; d++) {
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s).rows(0, M - 1) = softmax(
            gamma_B(d).slice(s) * X_B.slice(i).col(0)
          );
        }
        B(d).each_slice() = Btmp.t();
        log_B(d) = arma::log(B(d));
      }
    }
  }
  void update_probs(const arma::uword i) {
    update_pi(i);
    update_A(i);
    update_B(i);
  }
  void update_log_py(const arma::uword i) {
    for (arma::uword d = 0; d < D; d++) {
      for (arma::uword t = 0; t < Ti(i); t++) {
        log_py.slice(d).col(t) = log_B(d).slice(t).col(obs(t, i));
      }
    }
  }
};
#endif
