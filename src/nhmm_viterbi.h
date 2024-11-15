//Viterbi algorithm for NHMM and MHMM, single_sequence
#ifndef VITERBI_NHMM_H
#define VITERBI_NHMM_H

#include <RcppArmadillo.h>

template<typename subcol>
double univariate_viterbi_nhmm(
    subcol& q,
    const arma::vec& log_pi, 
    const arma::cube& log_A, 
    const arma::mat& log_py) {
  
  arma::uword S = log_py.n_rows;
  arma::uword T = log_py.n_cols;
  
  arma::mat delta(S, T);
  arma::umat phi(S, T);
  delta.col(0) = log_pi + log_py.col(0);
  phi.col(0).zeros();
  for (arma::uword t = 1; t < T; t++) {
    for (arma::uword j = 0; j < S; j++) {
      phi(j, t) = (delta.col(t - 1) + log_A.slice(t).col(j)).index_max();
      delta(j, t) = delta(phi(j, t), t - 1) + log_A(phi(j, t), j, t) + log_py(j, t);
    }
  }
  q(T - 1) = delta.col(T - 1).index_max();
  for (int t = (T - 2); t >= 0; t--) {
    q(t) = phi(q(t + 1), t + 1);
  }
  return delta.col(T - 1).max();
}

template <typename Model>
void viterbi_nhmm(Model& model, arma::umat& q, arma::vec& logp) {
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
    arma::subview_col<unsigned int> subcol = q.col(i);
    logp(i) = univariate_viterbi_nhmm(
      subcol,
      model.log_Pi, 
      model.log_A.slices(0, model.Ti(i) - 1), 
      model.log_py.cols(0, model.Ti(i) - 1)
    );
  }
}
template <typename Model>
void viterbi_mnhmm(Model& model, arma::umat& q, arma::vec& logp) {
  logp.fill(-arma::datum::inf);
  double logp_d;
  arma::uvec q_d(q.n_rows);
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
      logp_d = univariate_viterbi_nhmm(
        q_d,
        model.log_omega(d) + model.log_Pi(d),
        model.log_A(d).slices(0, model.Ti(i) - 1),
        model.log_py.slice(d).cols(0, model.Ti(i) - 1)
      );
      if (logp_d > logp(i)) {
        logp(i) = logp_d;
        q.col(i) = q_d;
      }
    }
  }
}

#endif
