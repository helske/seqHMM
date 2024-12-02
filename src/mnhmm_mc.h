#ifndef MNHMM_MC_H
#define MNHMM_MC_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "mnhmm_base.h"

struct mnhmm_mc : public mnhmm_base {
  
  const arma::ucube& obs;
  const arma::uword C;
  arma::field<arma::cube> eta_B;
  arma::uvec M;
  arma::field<arma::mat> Qm;
  arma::field<arma::cube> gamma_B;
  // these store pi, A, B, and log_p(y) of _one_ id and cluster we are currently working with
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  // excepted counts for EM algorithm
  arma::uword current_c; // for EM
  arma::uword current_d; // for EM
  arma::field<arma::cube> E_B;
  
  mnhmm_mc(
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
    const arma::ucube& obs_,
    const arma::mat& eta_omega_,
    const arma::field<arma::mat>& eta_pi_,
    const arma::field<arma::cube>& eta_A_,
    const arma::field<arma::cube>& eta_B_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0)
    : mnhmm_base(S_, D_, X_d_, X_pi_, X_s_, X_o_, Ti_, icpt_only_omega_, 
      icpt_only_pi_, icpt_only_A_, icpt_only_B_, iv_A_, 
      iv_B_, tv_A_, tv_B_, eta_omega_, eta_pi_, eta_A_, n_obs_, lambda_),
      obs(obs_),
      C(obs.n_rows),
      eta_B(C,D),
      M(C),
      Qm(C),
      gamma_B(C,D),
      B(C,D),
      log_B(C,D),
      current_c(0),
      current_d(0),
      E_B(C, D) {
    
    for (arma::uword d = 0; d < D; d++) {
      for (arma::uword c = 0; c < C; c++) {
        eta_B(c, d) = eta_B_(c + d * C); // 1D field from R to 2D field
        M(c) = eta_B(c, d).n_rows + 1;
        Qm(c) = create_Q(M(c));
        gamma_B(c, d) = eta_to_gamma(eta_B(c, d), Qm(c));
        B(c, d) = arma::cube(S, M(c) + 1, T);   // B field initialization
        log_B(c, d) = arma::cube(S, M(c) + 1, T); // log_B field initialization
        E_B(c, d) = arma::cube(T, N, S);
      }
    }
  }
  void update_gamma_B() {
    for (arma::uword c = 0; c < C; c++) {
      for (arma::uword d = 0; d < D; d++) {
        gamma_B(c, d) = eta_to_gamma(eta_B(c, d), Qm(c));
      }
    }
  }
  void update_B(const arma::uword i) {
    if (icpt_only_B) {
      for (arma::uword c = 0; c < C; c++) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        for (arma::uword d = 0; d < D; d++) {
          for (arma::uword s = 0; s < S; s++) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(c, d).slice(s).col(0)
            );
          }
          B(c, d).each_slice() = Btmp.t();
          log_B(c, d) = arma::log(B(c, d));
        }
      }
    } else { 
      if (tv_B) {
        for (arma::uword c = 0; c < C; c++) {
          arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
          for (arma::uword d = 0; d < D; d++) {
            for (arma::uword t = 0; t < Ti(i); t++) { // time
              for (arma::uword s = 0; s < S; s++) { // from states
                Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(c, d).slice(s) * X_B.slice(i).col(t));
              }
              B(c, d).slice(t) = Btmp.t();
            }
            log_B(c, d) = arma::log(B(c, d));
          }
        }
      } else {
        for (arma::uword c = 0; c < C; c++) {
          arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
          for (arma::uword d = 0; d < D; d++) {
            for (arma::uword s = 0; s < S; s++) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(
                gamma_B(c, d).slice(s) * X_B.slice(i).col(0)
              );
            }
            B(c, d).each_slice() = Btmp.t();
            log_B(c, d) = arma::log(B(c, d));
          }
        }
      }
    }
  }
  
  void update_log_py(const arma::uword i) {
    log_py.zeros();
    for (arma::uword d = 0; d < D; d++) {
      for (arma::uword t = 0; t < Ti(i); t++) {
        for (arma::uword c = 0; c < C; c++) {
          log_py.slice(d).col(t) += log_B(c, d).slice(t).col(obs(c, t, i));
        }
      }
    }
  }
  void estep_B(const arma::uword i, const arma::uword d, 
               const arma::mat& log_alpha, const arma::mat& log_beta, 
               const double ll, const double pseudocount = 0) {
    for (arma::uword k = 0; k < S; k++) { // state
      for (arma::uword t = 0; t < Ti(i); t++) { // time
        double pp = exp(log_alpha(k, t) + log_beta(k, t) - ll);
        for (arma::uword c = 0; c < C; c++) { // channel
          if (obs(c, t, i) < M(c)) {
            E_B(c, d)(t, i, k) = pp + pseudocount;
          } else {
            E_B(c, d)(t, i, k) = 0.0;
          }
        }
      }
    }
  }
  void mstep_B(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const arma::uword print_level);
  
  double objective_B(const arma::vec& x, arma::vec& grad);
  
  void compute_state_obs_probs(
      const arma::uword start, arma::field<arma::cube>& obs_prob, 
      arma::cube& state_prob
  );
};
#endif
