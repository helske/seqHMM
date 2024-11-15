#ifndef FORWARD_NHMM_H
#define FORWARD_NHMM_H

#include <RcppArmadillo.h>
#include "logsumexp.h"

// // time-varying A
template<typename submat>
void univariate_forward_nhmm(
    submat& log_alpha,
    const arma::vec& log_pi, 
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  log_alpha.col(0) = log_pi + log_py.col(0);
  for (arma::uword t = 1; t < T; t++) {
    for (arma::uword i = 0; i < S; i++) {
      log_alpha(i, t) = logSumExp(
        log_alpha.col(t - 1) + log_A.slice(t - 1).col(i) + log_py(i, t)
      );
    }
  }
}

template <typename Model>
void forward_nhmm(Model& model, arma::cube& log_alpha) {
  for (arma::uword i = 0; i < model.N; i++) {
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    univariate_forward_nhmm(
      log_alpha.slice(i),
      model.log_Pi,
      model.log_A, 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
  }
}

template <typename Model>
void forward_mnhmm(Model& model, arma::cube& log_alpha) {
  for (arma::uword i = 0; i < model.N; i++) {
    if (!model.icpt_only_omega || i == 0) {
      model.update_omega(i);
    }
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    for (arma::uword d = 0; d < model.D; d++) {
      arma::subview<double> submat = 
        log_alpha.slice(i).rows(d * model.S, (d + 1) * model.S - 1);
      univariate_forward_nhmm(
        submat,
        model.log_omega(d) + model.log_Pi(d),
        model.log_A(d), 
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
    }
  }
}
#endif
