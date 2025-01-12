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
  // these store pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  // excepted counts for EM algorithm
  arma::uword current_c; // for EM
  arma::field<arma::cube> E_B;
  
  nhmm_mc(
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
    const arma::ucube& obs_,
    const arma::mat& eta_pi_,
    const arma::cube& eta_A_,
    const arma::field<arma::cube>& eta_B_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0)
    : nhmm_base(S_, X_pi_, X_s_, X_o_, Ti_, icpt_only_pi_, icpt_only_A_, 
      icpt_only_B_, iv_A_, iv_B_, tv_A_, tv_B_, eta_pi_, eta_A_, n_obs_, 
      lambda_),
      obs(obs_), 
      C(obs.n_rows), 
      eta_B(eta_B_),
      M(C), 
      Qm(C),
      gamma_B(C),
      B(C),
      log_B(C),
      current_c(0),
      E_B(C) {
    for (arma::uword c = 0; c < C; c++) {
      M(c) = eta_B(c).n_rows + 1;
      Qm(c) = create_Q(M(c));
      gamma_B(c) = eta_to_gamma(eta_B(c), Qm(c));
      B(c) = arma::cube(S, M(c) + 1, T);   // B field initialization
      log_B(c) = arma::cube(S, M(c) + 1, T); // log_B field initialization
      E_B(c) = arma::cube(T, N, S);
    }
  }
  
  void update_gamma_B() {
    for (arma::uword c = 0; c < C; c++) {
      gamma_B(c) = eta_to_gamma(eta_B(c), Qm(c));
    }
  }
  
  void update_B(const arma::uword i) {
    if (icpt_only_B) {
      for (arma::uword c = 0; c < C; c++) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(c).slice(s).col(0)
          );
        }
        B(c).each_slice() = Btmp.t();
        log_B(c) = arma::log(B(c));
      }
    } else {
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
  }
  
  void update_log_py(const arma::uword i) {
    log_py.zeros();
    for (arma::uword t = 0; t < Ti(i); t++) {
      for (arma::uword c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
  }
  void estep_B(const arma::uword i, const arma::mat& log_alpha, 
               const arma::mat& log_beta, const double ll) {
    for (arma::uword k = 0; k < S; k++) { // state
      for (arma::uword t = 0; t < Ti(i); t++) { // time
        double pp = exp(log_alpha(k, t) + log_beta(k, t) - ll);
        for (arma::uword c = 0; c < C; c++) { // channel
          if (obs(c, t, i) < M(c) && pp > minval) {
            E_B(c)(t, i, k) = pp;
          } else {
            E_B(c)(t, i, k) = 0.0;
          }
        }
      }
    }
  }
  
  void mstep_B(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const double bound, 
               const arma::uword print_level);
  
  double objective_B(const arma::vec& x, arma::vec& grad);
  
  void compute_state_obs_probs(
      const arma::uword start, arma::field<arma::cube>& obs_prob, 
      arma::cube& state_prob
  );
};
#endif
