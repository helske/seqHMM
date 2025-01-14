#ifndef NHMM_SC_H
#define NHMM_SC_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "nhmm_base.h"

struct nhmm_sc : public nhmm_base {
  const arma::umat& obs;
  arma::cube eta_B;
  const arma::uword M;
  arma::mat Qm;
  arma::cube gamma_B;
  // these store pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::cube B;
  arma::cube log_B;
  // excepted counts for EM algorithm
  arma::cube E_B;
  
  nhmm_sc(
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
    const arma::umat& obs_,
    const arma::mat& eta_pi_,
    const arma::cube& eta_A_,
    const arma::cube& eta_B_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0)
    : nhmm_base(S_, X_pi_, X_s_, X_o_, Ti_, icpt_only_pi_, icpt_only_A_, 
      icpt_only_B_, iv_A_, iv_B_, tv_A_, tv_B_, eta_pi_, eta_A_, n_obs_, 
      lambda_),
      obs(obs_),  
      eta_B(eta_B_), 
      M(eta_B.n_rows + 1), 
      Qm(create_Q(M)), 
      gamma_B(eta_to_gamma(eta_B, Qm)),
      B(S, M + 1, T),
      log_B(S, M + 1, T), 
      E_B(T, N, S) {
  }
  
  void update_gamma_B() {
    gamma_B = eta_to_gamma(eta_B, Qm);
  }
  
  void update_B(const arma::uword i) {
    arma::mat Btmp(M + 1, S, arma::fill::ones);
    if (icpt_only_B) {
      for (arma::uword s = 0; s < S; s++) { // from states
        Btmp.col(s).rows(0, M - 1) = softmax(
          gamma_B.slice(s).col(0)
        );
      }
      B.each_slice() = Btmp.t();
    } else {
      if (tv_B) {
        for (arma::uword t = 0; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s++) { // from states
            Btmp.col(s).rows(0, M - 1) = softmax(gamma_B.slice(s) * X_B.slice(i).col(t));
          }
          B.slice(t) = Btmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s).rows(0, M - 1) = softmax(
            gamma_B.slice(s) * X_B.slice(i).col(0)
          );
        }
        B.each_slice() = Btmp.t();
      }
    }
    log_B = arma::log(B);
  }
  
  void update_log_py(const arma::uword i) {
    for (arma::uword t = 0; t < Ti(i); t++) {
      log_py.col(t) = log_B.slice(t).col(obs(t, i));
    }
  }
  void estep_B(const arma::uword i, const arma::mat& log_alpha, 
               const arma::mat& log_beta, const double ll) {
    for (arma::uword k = 0; k < S; k++) { // state
      for (arma::uword t = 0; t < Ti(i); t++) { // time
        if (obs(t, i) < M) {
          E_B(t, i, k) = exp(log_alpha(k, t) + log_beta(k, t) - ll);
        } else {
          E_B(t, i, k) = 0.0;
        }
      }
    }
    E_B.col(i).clean(minval);
  }
  
  void mstep_B(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const double bound, 
               const arma::uword print_level);
  
  double objective_B(const arma::vec& x, arma::vec& grad);
  
  void compute_state_obs_probs(
      const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob
  );
  void compute_state_obs_probs_fanhmm(
      const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob,
      const arma::field<arma::cube>& W_A, const arma::field<arma::cube>& W_B); 
};

#endif
