// simulate NHMMs
#include "get_parameters.h"
#include "eta_to_gamma.h"
#include "create_Q.h"
#include "sample.h"

// [[Rcpp::export]]
Rcpp::List simulate_nhmm_cpp(
    const arma::mat& eta_pi, const arma::mat& X_pi,
    const arma::cube& eta_A, const arma::cube& X_A,
    const arma::field<arma::cube>& eta_B, const arma::field<arma::cube>& X_B, 
    const arma::uvec& M) {
  arma::uword N = X_A.n_slices;
  arma::uword T = X_A.n_cols;
  arma::uword S = eta_A.n_slices;
  arma::uword C = M.n_elem;
  
  arma::ucube y(C, T, N);
  arma::umat z(T, N);
  arma::mat Qs = create_Q(S);
  arma::field<arma::mat> Qm = create_Q(M);
  arma::mat gamma_pi = eta_to_gamma(eta_pi, Qs);
  arma::cube gamma_A = eta_to_gamma(eta_A, Qs);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B, Qm);
  arma::vec pi(S);
  arma::cube A(S, S, T);
  arma::field<arma::cube> B(C);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::field<arma::uvec> seqM(C);
  for (arma::uword c = 0; c < C; ++c) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  for (arma::uword i = 0; i < N; ++i) {
    pi = get_pi(gamma_pi, X_pi.col(i));
    A = get_A(gamma_A, X_A.slice(i));
    z(0, i) = arma::as_scalar(sample(seqS, pi));
    for (arma::uword c = 0; c < C; ++c) {
      B(c) = get_B(gamma_B(c), X_B(c).slice(i), M(c));
      y(c, 0, i) = arma::as_scalar(sample(seqM(c), B(c).slice(0).row(z(0, i)).t()));
    }
    for (arma::uword t = 1; t < T; ++t) {
      z(t, i) = arma::as_scalar(sample(seqS, A.slice(t).row(z(t - 1, i)).t()));
      for (arma::uword c = 0; c < C; ++c) {
        y(c, t, i) = arma::as_scalar(sample(seqM(c), B(c).slice(t).row(z(t, i)).t()));
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y), 
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}


// // [[Rcpp::export]]
// Rcpp::List simulate_fanhmm_cpp(
//     const arma::mat& eta_pi, const arma::mat& X_pi,
//     const arma::cube& eta_A, const arma::field<arma::cube>& X_A,
//     const Rcpp::List& eta_B, const Rcpp::List& X_B_,
//     const arma::umat& obs_1, const bool autoregression, 
//     const arma::uvec& M) {
//   arma::uword N = obs_1.n_elem;
//   arma::uword T = X_A(0).n_cols;
//   arma::uword S = eta_A.n_slices;
//   arma::uword C = M.n_elem;
//   arma::ucube y(C, T, N);
//   arma::umat z(T, N);
//   arma::mat Qs = create_Q(S);
//   arma::field<arma::mat> Qm = create_Q(M);
//   arma::mat gamma_pi = eta_to_gamma(eta_pi, Qs);
//   arma::cube gamma_A = eta_to_gamma(eta_A, Qs);
//   arma::field<arma::cube> gamma_B = eta_to_gamma(list_to_2d_field(eta_B), Qm);
//   arma::vec pi(S);
//   arma::vec A(S);
//   arma::field<arma::vec> B(C);
//   arma::field<arma::uvec> seqM(C);
//   for (arma::uword c = 0; c < C; ++c) {
//     seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
//   }
//   arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
//   arma::field<arma::cube> X_B = list_to_2d_field(X_B_); // C x M
//   for (arma::uword i = 0; i < N; ++i) {
//     pi = get_pi(gamma_pi, X_pi.col(i));
//     z(0, i) = arma::as_scalar(sample(seqS, pi));
//     if (autoregression) {
//       for (arma::uword c = 0; c < C; ++c) {
//         y(c, 0, i) = obs_1(i, c);
//       }
//     } else {
//       for (arma::uword c = 0; c < C; ++c) {
//         B(c) = softmax(
//           gamma_B(c).slice(z(0, i)) * X_B(0, c).slice(i).col(0)
//         );
//         y(c, 0, i) = arma::as_scalar(sample(seqM(c), B(c)));
//       }
//     }
//     for (arma::uword t = 1; t < T; ++t) {
//       for (arma::uword c = 0; c < C; ++c) {
//         A = softmax(
//           gamma_A.slice(z(t - 1, i)) * X_A(y(c, t - 1, i)).slice(i).col(t)
//         );
//         z(t, i) = arma::as_scalar(sample(seqS, A));
//         B(c) = softmax(
//           gamma_B(c).slice(z(t, i)) * X_B(y(c, t - 1, i)).slice(i).col(t)
//         );
//         y(c, t, i) = arma::as_scalar(sample(seqM(c), B(c)));
//       }
//     }
//   }
//     return Rcpp::List::create(
//       Rcpp::Named("observations") = Rcpp::wrap(y), 
//       Rcpp::Named("states") = Rcpp::wrap(z)
//     );
//   }
  