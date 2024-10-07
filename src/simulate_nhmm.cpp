// simulate NHMMs
#include <RcppArmadilloExtensions/sample.h>
#include "get_parameters.h"
#include "eta_to_gamma.h"

// [[Rcpp::export]]
Rcpp::List simulate_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::cube& eta_B, const arma::cube& X_o) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  arma::umat y(T, N);
  arma::umat z(T, N);
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::cube gamma_B = eta_to_gamma(eta_B);
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::cube B(S, M, T);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqM = arma::linspace<arma::uvec>(0, M - 1, M);
  for (unsigned int i = 0; i < N; i++) {
    Pi = get_pi(gamma_pi, X_i.col(i));
    A = get_A(gamma_A, X_s.slice(i));
    B = get_B(gamma_B, X_o.slice(i), false);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, Pi));
    y(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(0).row(z(0, i)).t()));
    for (unsigned int t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      y(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(t).row(z(t, i)).t()));
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
// [[Rcpp::export]]
Rcpp::List simulate_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o, 
    const arma::uvec& M) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int C = M.n_elem;
  
  arma::ucube y(C, T, N);
  arma::umat z(T, N);

  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::field<arma::uvec> seqM(C);
  for (unsigned int c = 0; c < C; c++) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
    B(c) = arma::cube(M(c) - 1, X_o.n_rows, S);
  }
  for (unsigned int i = 0; i < N; i++) {
    Pi = get_pi(gamma_pi, X_i.col(i));
    A = get_A(gamma_A, X_s.slice(i));
    B = get_B(gamma_B, X_o.slice(i), M, false);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, Pi));
    for (unsigned int c = 0; c < C; c++) {
      y(c, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(0).row(z(0, i)).t()));
    }
    for (unsigned int t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      for (unsigned int c = 0; c < C; c++) {
        y(c, t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(t).row(z(t, i)).t()));
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}

// [[Rcpp::export]]
Rcpp::List simulate_mnhmm_singlechannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  unsigned int D = eta_omega.n_rows + 1;
  arma::umat y(T, N);
  arma::umat z(T, N);
  
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec omega(D);
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::cube B(S, M, T);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqM = arma::linspace<arma::uvec>(0, M - 1, M);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  for (unsigned int i = 0; i < N; i++) {
    omega = get_omega(gamma_omega, X_d.col(i));
    unsigned int cluster = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqD, 1, false, omega));
    Pi = get_pi(gamma_pi(cluster), X_i.col(i));
    A = get_A(gamma_A(cluster), X_s.slice(i));
    B = get_B(gamma_B(cluster), X_o.slice(i), false);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, Pi));
    y(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(0).row(z(0, i)).t()));
    for (unsigned int t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      y(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(t).row(z(t, i)).t()));
    }
    z.col(i) += cluster * S;
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}


// [[Rcpp::export]]
Rcpp::List simulate_mnhmm_multichannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d, 
    const arma::uvec& M) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int C = M.n_elem;
  arma::ucube y(C, T, N);
  arma::umat z(T, N);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec omega(D);
  arma::vec Pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::field<arma::uvec> seqM(C);
  for (unsigned int c = 0; c < C; c++) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
    B(c) = arma::cube(M(c) - 1, X_o.n_rows, S);
  }
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  for (unsigned int i = 0; i < N; i++) {
    omega = get_omega(gamma_omega, X_d.col(i));
    unsigned int cluster = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqD, 1, false, omega));
    Pi = get_pi(gamma_pi(cluster), X_i.col(i));
    A = get_A(gamma_A(cluster), X_s.slice(i));
    B = get_B(gamma_B.rows(cluster * C, (cluster + 1) * C - 1), X_o.slice(i), M, false);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, Pi));
    for (unsigned int c = 0; c < C; c++) {
      y(c, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(0).row(z(0, i)).t()));
    }
    for (unsigned int t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      for (unsigned int c = 0; c < C; c++) {
        y(c, t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(t).row(z(t, i)).t()));
      }
    }
    z.col(i) += cluster * S;
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}