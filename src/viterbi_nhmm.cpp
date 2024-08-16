//Viterbi algorithm for NHMM and MHMM, single_sequence
#include "get_parameters.h"
#include "viterbi_nhmm.h"

// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_singlechannel(const arma::mat& beta_i_raw, const arma::mat& X_i,
                        const arma::cube& beta_s_raw, const arma::cube& X_s,
                        const arma::cube& beta_o_raw, const arma::cube& X_o,
                        const arma::mat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  arma::umat q(T, N);
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
Rcpp::List viterbi_nhmm_multichannel(const arma::mat& beta_i_raw, const arma::mat& X_i,
                                      const arma::cube& beta_s_raw, const arma::cube& X_s,
                                      const arma::vec& beta_o_raw, const arma::cube& X_o,
                                      const arma::cube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = beta_s_raw.n_slices;
  unsigned int C = obs.n_rows;
  arma::umat q(T, N);
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

double univariate_viterbi_nhmm(const arma::vec& init, const arma::cube& transition, 
                               const arma::mat& log_py, arma::subview_col<unsigned int> q) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  arma::mat delta(S, T);
  arma::umat phi(S, T);
  delta.col(0) = init + log_py.col(0);
  phi.col(0).zeros();
  for (unsigned int t = 1; t < T; t++) {
    for (unsigned int j = 0; j < S; j++) {
      phi(j, t) = (delta.col(t - 1) + transition.slice(t).col(j)).index_max();
      delta(j, t) = delta(phi(j, t), t - 1) + transition(phi(j, t), j, t) + log_py(j, t);
    }
  }
  q(T - 1) = delta.col(T - 1).index_max();
  for (int t = (T - 2); t >= 0; t--) {
    q(t) = phi(q(t + 1), t + 1);
  }
  return delta.col(T - 1).max();
}
