#ifndef FANHMM_H
#define FANHMM_H

#include "config.h"
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "nhmm.h"

struct fanhmm : public nhmm {
  
  const arma::field<arma::cube>& W;
  const arma::field<arma::vec>& delta;
  fanhmm(
    const arma::ucube& obs_,
    const arma::uvec& Ti_,
    const arma::uvec& M_,
    const arma::mat& X_pi_,
    const arma::cube& X_A_,
    const arma::field<arma::cube>& X_B_,
    const bool icpt_only_pi_,
    const bool icpt_only_A_,
    const arma::uvec& icpt_only_B_,
    const bool iv_A_,
    const arma::uvec& iv_B_,
    const bool tv_A_,
    const arma::uvec& tv_B_,
    const arma::mat& eta_pi_,
    const arma::cube& eta_A_,
    const arma::field<arma::cube>& eta_B_,
    const arma::field<arma::cube>& W_,
    const arma::field<arma::vec>& delta_,
    const double lambda_ = 0,
    double maxval_ = arma::datum::inf,
    double minval_ = -1.0)
    : nhmm(obs_, Ti_, M_, X_pi_, X_A_, X_B_, icpt_only_pi_, icpt_only_A_, 
      icpt_only_B_, iv_A_, iv_B_, tv_A_, tv_B_, eta_pi_, eta_A_, eta_B_, 
      lambda_, maxval_, minval_), W(W_), delta(delta_) {
  }
  
  void update_B(const arma::uword i) {
    for (arma::uword c = 0; c < C; ++c) {
      if (icpt_only_B(c)) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(c).slice(s).col(0)
          );
        }
        B(c).each_slice() = Btmp.t();
        log_B(c) = arma::log(B(c));
      } else {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        if (tv_B(c)) {
          Btmp.cols(0, M(c) - 1).zeros();
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(c).slice(s) * W(c).slice(i).col(0)
            ) * delta(0);
          }
          for (arma::uword m = 1; m < M(c); ++m) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) += softmax(
                gamma_B(c).slice(s) * W(c).slice(i).col(m)
              ) * delta(c, m);
            }
          }
          B(c).slice(0) = Btmp.t();
          for (arma::uword t = 1; t < Ti(i); ++t) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(c).slice(s) * X_B(c).slice(i).col(t));
            }
            B(c).slice(t) = Btmp.t();
          }
        } else {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(c).slice(s) * X_B(c).slice(i).col(0)
            );
          }
          B(c).each_slice() = Btmp.t();
        }
        log_B(c) = arma::log(B(c));
      }
    }
  }
  void predict(
      arma::field<arma::cube>& obs_prob, const arma::field<arma::cube>& W_A, 
      const arma::field<arma::cube>& W_B
  );
};

#endif
