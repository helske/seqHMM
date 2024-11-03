#ifndef NHMM_MC_H
#define NHMM_MC_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "nhmm_base.h"

struct nhmm_mc : public nhmm_base {
  
  const arma::ucube& obs;
  const arma::uword C;
  arma::field<arma::cube> eta_B;
  arma::uvec M;
  arma::field<arma::mat> Qm;
  arma::field<arma::cube> gamma_B;
  // these store Pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  
  nhmm_mc(
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
    const arma::ucube& obs_,
    arma::mat& eta_pi_,
    arma::cube& eta_A_,
    arma::field<arma::cube>& eta_B_)
    : nhmm_base(S_, X_pi_, X_s_, X_o_, Ti_, iv_pi_, iv_A_, iv_B_, tv_A_, tv_B_, eta_pi_, eta_A_),
      obs(obs_), 
      C(obs.n_rows), 
      eta_B(eta_B_),
      M(arma::uvec(C)), 
      Qm(arma::field<arma::mat>(C)),
      gamma_B(arma::field<arma::cube>(C)),
      B(arma::field<arma::cube>(C)),
      log_B(arma::field<arma::cube>(C)) {
    
    for (arma::uword c = 0; c < C; c++) {
      M(c) = eta_B(c).n_rows + 1;
      Qm(c) = create_Q(M(c));
      gamma_B(c) = eta_to_gamma(eta_B(c), Qm(c));
      B(c) = arma::cube(S, M(c) + 1, T);   // B field initialization
      log_B(c) = arma::cube(S, M(c) + 1, T); // log_B field initialization
    }
  }
  void update_B(const arma::uword i) {
    for (arma::uword c = 0; c < C; c++) {
      arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
      if (tv_B) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s++) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(c).slice(s) * X_B.slice(i).col(t));
          }
          B(c).slice(t) = Btmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(c).slice(s) * X_B.slice(i).col(0)
          );
        }
        B(c).each_slice() = Btmp.t();
      }
      log_B(c) = arma::log(B(c));
    }
  }
  
  void update_probs(const arma::uword i) {
    update_pi(i);
    update_A(i);
    update_B(i);
  }
  
  void update_log_py(const arma::uword i) {
    log_py.zeros();
    for (arma::uword t = 0; t < Ti(i); t++) {
      for (arma::uword c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
  }
};
#endif
