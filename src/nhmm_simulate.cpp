// simulate NHMMs
#include <RcppArmadilloExtensions/sample.h>
#include "get_parameters.h"
#include "eta_to_gamma.h"

#include "create_Q.h"

// [[Rcpp::export]]
Rcpp::List simulate_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_pi,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::cube& eta_B, const arma::cube& X_B) {
  arma::uword N = X_A.n_slices;
  arma::uword T = X_A.n_cols;
  arma::uword S = eta_A.n_slices;
  arma::uword M = eta_B.n_rows + 1;
  arma::ucube y(1, T, N);
  arma::umat z(T, N);
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::cube gamma_B = eta_to_gamma(eta_B);
  arma::vec pi(S);
  arma::cube A(S, S, T);
  arma::cube B(S, M, T);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqM = arma::linspace<arma::uvec>(0, M - 1, M);
  for (arma::uword i = 0; i < N; i++) {
    pi = get_pi(gamma_pi, X_pi.col(i));
    A = get_A(gamma_A, X_A.slice(i));
    B = get_B(gamma_B, X_B.slice(i));
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, pi));
    y(0, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(0).row(z(0, i)).t()));
    for (arma::uword t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      y(0, t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(t).row(z(t, i)).t()));
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
// [[Rcpp::export]]
Rcpp::List simulate_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& X_pi,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_B, 
    const arma::uvec& M) {
  arma::uword N = X_A.n_slices;
  arma::uword T = X_A.n_cols;
  arma::uword S = eta_A.n_slices;
  arma::uword C = M.n_elem;
  
  arma::ucube y(C, T, N);
  arma::umat z(T, N);
  
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::field<arma::uvec> seqM(C);
  for (arma::uword c = 0; c < C; c++) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  for (arma::uword i = 0; i < N; i++) {
    pi = get_pi(gamma_pi, X_pi.col(i));
    A = get_A(gamma_A, X_A.slice(i));
    B = get_B(gamma_B, X_B.slice(i), M);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, pi));
    for (arma::uword c = 0; c < C; c++) {
      y(c, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(0).row(z(0, i)).t()));
    }
    for (arma::uword t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      for (arma::uword c = 0; c < C; c++) {
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
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::mat& eta_omega, const arma::mat& X_omega) {
  arma::uword N = X_A.n_slices;
  arma::uword T = X_A.n_cols;
  arma::uword S = eta_A(0).n_slices;
  arma::uword M = eta_B(0).n_rows + 1;
  arma::uword D = eta_omega.n_rows + 1;
  arma::ucube y(1, T, N);
  arma::umat z(T, N);
  
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec omega(D);
  arma::vec pi(S);
  arma::cube A(S, S, T);
  arma::cube B(S, M, T);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqM = arma::linspace<arma::uvec>(0, M - 1, M);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  for (arma::uword i = 0; i < N; i++) {
    omega = get_omega(gamma_omega, X_omega.col(i));
    arma::uword cluster = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqD, 1, false, omega));
    pi = get_pi(gamma_pi(cluster), X_pi.col(i));
    A = get_A(gamma_A(cluster), X_A.slice(i));
    B = get_B(gamma_B(cluster), X_B.slice(i));
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, pi));
    y(0, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(0).row(z(0, i)).t()));
    for (arma::uword t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      y(0, t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B.slice(t).row(z(t, i)).t()));
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
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_pi,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_A,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_B,
    const arma::mat& eta_omega, const arma::mat& X_omega, 
    const arma::uvec& M) {
  arma::uword N = X_A.n_slices;
  arma::uword T = X_A.n_cols;
  arma::uword S = eta_A(0).n_slices;
  arma::uword D = eta_omega.n_rows + 1;
  arma::uword C = M.n_elem;
  arma::ucube y(C, T, N);
  arma::umat z(T, N);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  arma::vec omega(D);
  arma::vec pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::field<arma::uvec> seqM(C);
  for (arma::uword c = 0; c < C; c++) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  for (arma::uword i = 0; i < N; i++) {
    omega = get_omega(gamma_omega, X_omega.col(i));
    arma::uword cluster = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqD, 1, false, omega));
    pi = get_pi(gamma_pi(cluster), X_pi.col(i));
    A = get_A(gamma_A(cluster), X_A.slice(i));
    B = get_B(gamma_B.rows(cluster * C, (cluster + 1) * C - 1), X_B.slice(i), M);
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, pi));
    for (arma::uword c = 0; c < C; c++) {
      y(c, 0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM(c), 1, false, B(c).slice(0).row(z(0, i)).t()));
    }
    for (arma::uword t = 1; t < T; t++) {
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A.slice(t).row(z(t - 1, i)).t()));
      for (arma::uword c = 0; c < C; c++) {
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

// [[Rcpp::export]]
Rcpp::List simulate_fanhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_pi,
    const arma::cube& eta_A, const arma::field<arma::cube>& X_A,
    const arma::cube& eta_B, const arma::field<arma::cube>& X_B,
    const arma::uvec& obs_1, const bool autoregression) {
  arma::uword N = obs_1.n_elem;
  arma::uword T = X_A(0).n_cols;
  arma::uword S = eta_A.n_slices;
  arma::uword M = eta_B.n_rows + 1;
  arma::umat y(T, N);
  arma::umat z(T, N);
  arma::mat Qs = create_Q(S);
  arma::mat Qm = create_Q(M);
  arma::mat gamma_pi = eta_to_gamma(eta_pi, Qs);
  arma::cube gamma_A = eta_to_gamma(eta_A, Qs);
  arma::cube gamma_B = eta_to_gamma(eta_B, Qm);
  arma::vec pi(S);
  arma::vec A(S);
  arma::vec B(S);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqM = arma::linspace<arma::uvec>(0, M - 1, M);
  
  for (arma::uword i = 0; i < N; i++) {
    pi = get_pi(gamma_pi, X_pi.col(i));
    z(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, pi));
   
    if (autoregression) {
      y(0, i) = obs_1(i);
    } else {
      B = softmax(
        gamma_B.slice(z(0, i)) * X_B(0).slice(i).col(0)
      );
      y(0, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B));
    }
    for (arma::uword t = 1; t < T; t++) {
      A = softmax(
        gamma_A.slice(z(t - 1, i)) * X_A(y(t - 1, i)).slice(i).col(t)
      );
      z(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqS, 1, false, A));
      B = softmax(
        gamma_B.slice(z(t, i)) * X_B(y(t - 1, i)).slice(i).col(t)
      );
      y(t, i) = arma::as_scalar(Rcpp::RcppArmadillo::sample(seqM, 1, false, B));
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
