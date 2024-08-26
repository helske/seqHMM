//Viterbi algorithm for NHMM and MHMM, single_sequence
#include "get_parameters.h"
#include "viterbi_nhmm.h"

double univariate_viterbi_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py,
    arma::subview_col<unsigned int> q) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  arma::mat delta(S, T);
  arma::umat phi(S, T);
  delta.col(0) = log_init + log_py.col(0);
  phi.col(0).zeros();
  for (unsigned int t = 1; t < T; t++) {
    for (unsigned int j = 0; j < S; j++) {
      phi(j, t) = (delta.col(t - 1) + log_transition.slice(t).col(j)).index_max();
      delta(j, t) = delta(phi(j, t), t - 1) + log_transition(phi(j, t), j, t) + log_py(j, t);
    }
  }
  q(T - 1) = delta.col(T - 1).index_max();
  for (int t = (T - 2); t >= 0; t--) {
    q(t) = phi(q(t + 1), t + 1);
  }
  return delta.col(T - 1).max();
}
// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_singlechannel(
    const arma::mat& beta_i_raw, const arma::mat& X_i,
    const arma::cube& beta_s_raw, const arma::cube& X_s,
    const arma::cube& beta_o_raw, const arma::cube& X_o,
    const arma::mat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(S, T);
  arma::mat log_Pi = get_pi(beta_i_raw, X_i, 1);
  // field of S x S x T cubes
  arma::field<arma::cube> log_A = get_A(beta_s_raw, X_s, 1);
  // field of S x M x T cubes
  arma::field<arma::cube> log_B = get_B(beta_o_raw, X_o, 1);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(i).at(s, obs(t, i), t);
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi.col(i), log_A(i), log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}

// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_multichannel(
    const arma::mat& beta_i_raw, const arma::mat& X_i,
    const arma::cube& beta_s_raw, const arma::cube& X_s,
    const arma::vec& beta_o_raw, const arma::cube& X_o,
    const arma::cube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  unsigned int C = obs.n_rows;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_Pi = get_pi(beta_i_raw, X_i, 1);
  arma::field<arma::cube> log_A = get_A(beta_s_raw, X_s, 1);
  arma::field<arma::cube> log_B = get_multichannel_B(beta_o_raw, X_o, S, C, M, 1);
  arma::mat log_py(S, T);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = 0;
        for (unsigned int c = 0; c < C; c++) {
          log_py(s, t) += log_B(i, c).at(s, obs(c, t, i), t);
        }
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi.col(i), log_A(i), log_py, q.col(i));
  }
  
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}

// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_singlechannel(
    const arma::mat& beta_i_raw, const arma::mat& X_i,
    const arma::cube& beta_s_raw, const arma::cube& X_s,
    const arma::cube& beta_o_raw, const arma::cube& X_o,
    const arma::mat& theta_raw, const arma::mat& X_d,
    const arma::mat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  unsigned int D = theta_raw.n_rows + 1;
  unsigned int SD = S * D;
  unsigned int M = beta_o_raw.n_rows + 1;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(SD, T);
  arma::vec log_Pi(SD);
  arma::cube log_A(SD, SD, T, arma::fill::value(-arma::datum::inf));
  arma::cube log_B(SD, M, T, arma::fill::value(-arma::datum::inf));
  for (unsigned int i = 0; i < N; i++) {
    arma::vec log_omega = get_omega_i(theta_raw, X_d.col(i), 1);
    for (unsigned int d = 0; d < D; d++) {
      log_Pi.rows(d * S, (d + 1) * S - 1) = log_omega(d) + 
        get_pi_i(beta_i_raw.rows(d * (S - 1), (d + 1) * (S - 1) - 1), X_i.col(i), 1);
      log_A.tube(d * S, d * S, (d + 1) * S - 1, (d + 1) * S - 1) = 
        get_A_i(beta_s_raw.rows(d * (S - 1), (d + 1) * (S - 1) - 1), X_s.slice(i), 1);
      log_B.rows(d * S, (d + 1) * S - 1) = 
        get_B_i(beta_o_raw.rows(d * (S - 1), (d + 1) * (S - 1) - 1), X_o.slice(i), 1);
    }
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(s, obs(t, i), t);
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_multichannel(
    const arma::mat& beta_i_raw, const arma::mat& X_i,
    const arma::cube& beta_s_raw, const arma::cube& X_s,
    const arma::mat& beta_o_raw, const arma::cube& X_o,
    const arma::mat& theta_raw, const arma::mat& X_d,
    const arma::cube& obs, const arma::uvec M) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  unsigned int D = theta_raw.n_rows + 1;
  unsigned int SD = S * D;
  unsigned int C = obs.n_rows;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(SD, T);
  arma::vec log_Pi(SD);
  arma::cube log_A(SD, SD, T, arma::fill::value(-arma::datum::inf));
  arma::field<arma::cube> log_B(C);
  for (unsigned int i = 0; i < N; i++) {
    arma::vec log_omega = get_omega_i(theta_raw, X_d.col(i), 1);
    for (unsigned int d = 0; d < D; d++) {
      log_Pi.rows(d * S, (d + 1) * S - 1) = log_omega(d) +
        get_pi_i(
          beta_i_raw.rows(d * (S - 1), (d + 1) * (S - 1) - 1), X_i.col(i), 1
        );
      log_A.tube(d * S, d * S, (d + 1) * S - 1, (d + 1) * S - 1) = get_A_i(
        beta_s_raw.rows(d * (S - 1), (d + 1) * (S - 1) - 1), X_s.slice(i), 1
      );
      log_B = get_multichannel_B_i(
        beta_o_raw.col(d), X_o.slice(i), S, C, M, 1
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(d * S + s, t) = 0;
          for (unsigned int c = 0; c < C; c++) {
            log_py(d * S + s, t) += log_B(c).at(s, obs(c, t, i), t);
          }
        }
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
