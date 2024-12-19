#ifndef FANHMM_SC_H
#define FANHMM_SC_H

#include <RcppArmadillo.h>
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "nhmm_sc.h"

struct fanhmm_sc : public nhmm_sc {
  
  const arma::uvec obs_0; // initial observation for each sequence
  const arma::cube& W_A;
  const arma::cube& W_B;
  const arma::uword L_A;
  const arma::uword L_B;
  // effects y_{t} on A
  // field of length S of (S-1) x L_A x (M-1) cubes (last class is reference)
  // rho_A(s) gives main and interaction effects of y_t when in state s
  arma::field<arma::cube> rho_A; 
  // effects y_{t-1} on B
  // field of length S of (M-1) x L_B x (M-1) cubes
  // rho_B(s) gives main and interaction effects of y_{t-1} when in state s
  arma::field<arma::cube> rho_B; 
  // effects y_{t} on A
  arma::field<arma::cube> phi_A;
  // effects y_{t-1} on B
  arma::field<arma::cube> phi_B;
  fanhmm_sc(
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
    const arma::uvec obs_0_,
    const arma::cube& W_A_,
    const arma::cube& W_B_,
    const arma::field<arma::cube>& rho_A_,
    const arma::field<arma::cube>& rho_B_,
    const arma::uword n_obs_ = 0,
    const double lambda_ = 0)
    : nhmm_sc(S_, X_pi_, X_s_, X_o_, Ti_, icpt_only_pi_, icpt_only_A_, 
      icpt_only_B_, iv_A_, iv_B_, tv_A_, tv_B_, obs_, eta_pi_, eta_A_, eta_B_, 
      n_obs_, lambda_),
      obs_0(obs_0_),
      W_A(W_A_),
      W_B(W_B_),
      L_A(W_A.n_rows),
      L_B(W_B.n_rows),
      rho_A(rho_A_),
      rho_B(rho_B_),
      phi_A(rho_to_phi(rho_A, Qs)),
      phi_B(rho_to_phi(rho_B, Qm)) {
  }
  void update_phi_A() {
    phi_A = rho_to_phi(rho_A, Qs);
  }
  void update_phi_B() {
    phi_B = rho_to_phi(rho_B, Qm);
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
            Atmp.col(s) = softmax(
              gamma_A.slice(s) * X_A.slice(i).col(t) +
                phi_A(s).slice(obs(t, i)) * W_A.slice(i).col(t)
            );
          }
          A.slice(t) = Atmp.t();
        }
      } else {
        // does not depend on time => does not depend on y_t
        for (arma::uword s = 0; s < S; s ++) { // from states
          Atmp.col(s) = softmax(gamma_A.slice(s) * X_A.slice(i).col(0));
        }
        A.each_slice() = Atmp.t();
      }
    }
    log_A = arma::log(A);
  }
  void update_B(const arma::uword i) {
    arma::mat Btmp(M, S);
    if (icpt_only_B) {
      for (arma::uword s = 0; s < S; s++) { // from states
        Btmp.col(s) = softmax(gamma_B.slice(s).col(0));
      }
      B.each_slice() = Btmp.t();
    } else {
      if (tv_B) {
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s) = 
            softmax(
              gamma_B.slice(s) * X_B.slice(i).col(0) +
                phi_B(s).slice(obs_0(i)) * W_B.slice(i).col(0)
            );
        }
        B.slice(0) = Btmp.t();
        for (arma::uword t = 1; t < Ti(i); t++) { // time
          for (arma::uword s = 0; s < S; s++) { // from states
            Btmp.col(s) = 
              softmax(
                gamma_B.slice(s) * X_B.slice(i).col(t) +
                  phi_B(s).slice(obs(t - 1, i)) * W_B.slice(i).col(t)
              );
          }
          B.slice(t) = Btmp.t();
        }
      } else {
        // does not depend on time => does not depend on y_t
        for (arma::uword s = 0; s < S; s++) { // from states
          Btmp.col(s) = softmax(
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
        E_B(t, i, k) = exp(log_alpha(k, t) + log_beta(k, t) - ll);
      }
    }
    E_B.col(i).clean(std::numeric_limits<double>::min());
  }
  
  void mstep_A(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const double bound, 
               const arma::uword print_level);
  
  void mstep_B(const double ftol_abs, const double ftol_rel, 
               const double xtol_abs, const double xtol_rel, 
               const arma::uword maxeval, const double bound, 
               const arma::uword print_level);
  
  double objective_A(const arma::vec& x, arma::vec& grad);
  double objective_B(const arma::vec& x, arma::vec& grad);
  
  void compute_state_obs_probs(
      const arma::uword start, arma::cube& obs_prob, arma::cube& state_prob
  );
};
#endif
