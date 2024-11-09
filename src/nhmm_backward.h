#ifndef BACKWARD_NHMM_H
#define BACKWARD_NHMM_H

#include <RcppArmadillo.h>
#include "logsumexp.h"

// time-varying A
template<typename submat>
void univariate_backward_nhmm(
    submat& log_beta,
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  
  log_beta.col(T - 1).zeros();
  for (int t = (T - 2); t >= 0; t--) {
    for (arma::uword i = 0; i < S; i++) {
      log_beta(i, t) = logSumExp(
        log_beta.col(t + 1) + log_A.slice(t).row(i).t() + 
          log_py.col(t + 1)
      );
    }
  }
}
// // time-invariant A
// template<typename submat>
// void univariate_backward_nhmm(
//     submat& log_beta,
//     const arma::mat& log_A, 
//     const arma::mat& log_py) {
//   
//   arma::uword S = log_py.n_rows;
//   arma::uword T = log_py.n_cols;
//   arma::mat log_tA = log_A.t();
//   log_beta.col(T - 1).zeros();
//   for (int t = (T - 2); t >= 0; t--) {
//     for (arma::uword i = 0; i < S; i++) {
//       log_beta(i, t) = logSumExp(
//         log_beta.col(t + 1) + log_tA.col(i) + log_py.col(t + 1)
//       );
//     }
//   }
// }

template <typename Model>
void backward_nhmm(Model& model, arma::cube& log_beta) {
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
    univariate_backward_nhmm(
      log_beta.slice(i),
      model.log_A, 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
  }
}


template <typename Model>
void backward_mnhmm(Model& model, arma::cube& log_beta) {
  for (arma::uword i = 0; i < model.N; i++) {
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (model.iv_B || i == 0) {
      model.update_B(i);
    }
    model.update_log_py(i);
    for (arma::uword d = 0; d < model.D; d++) {
      arma::subview<double> submat = 
        log_beta.slice(i).rows(d * model.S, (d + 1) * model.S - 1);
      univariate_backward_nhmm(
        submat,
        model.log_A(d), 
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
    }
  }
}

#endif
