// // EM algorithm for HMMs, used for initial fit of NHMMs
// 
// #include "eta_to_gamma.h"
// #include "softmax.h"
// #include "forward_nhmm.h"
// #include "backward_nhmm.h"
// 
// // [[Rcpp::export]]
// Rcpp::List EM_nhmm_singlechannel(
//     arma::mat& eta_pi, const arma::mat& X_pi,
//     arma::cube& eta_A, const arma::cube& X_A,
//     arma::cube& eta_B, const arma::cube& X_B,
//     const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
//     const bool tv_A, const bool tv_B, const arma::uvec& Ti,
//     arma::uword maxeval, double ftol_abs, double ftol_rel) {
//   
//   nhmm_sc model(
//       eta_A.n_slices, X_pi, X_A, X_B, Ti,
//       iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
//   );
//   
//   // EM-algorithm begins
//   
//   double relative_change = ftol_rel + 1.0;
//   double absolute_change = ftol_abs + 1.0;
//   arma::uword iter = 0;
//   double ll_new = 0;
//   double ll = -1e150;
//   // initial values to probabilities
//   arma::vec log_Pi = arma::log(softmax(eta_to_gamma(eta_pi)));
//   arma::mat log_A(S, S);
//   arma::mat log_B(S, M + 1);
//   log_B.col(M).zeros();
//   arma::cube gamma_A = eta_to_gamma(eta_A);
//   arma::cube gamma_B = eta_to_gamma(eta_B);
//   for (arma::uword s = 0; s < S; s++) {
//     log_A.row(s) = arma::log((softmax(gamma_A.slice(s)).t()));
//     log_B.row(s).cols(0, M - 1) = arma::log(softmax(gamma_B.slice(s)).t());
//   }
//   arma::vec Pi(S);
//   arma::mat A(S, S);
//   arma::mat B(S, M);
//   arma::mat log_alpha(S, T);
//   arma::mat log_beta(S, T);
//   arma::mat log_py(S, T);
//   arma::vec tmp_T(T);
//   while ((relative_change > ftol_rel) && (absolute_change > ftol_abs) && (iter < maxeval)) {
//     iter++;
//     ll_new = 0;
//     Pi.zeros();
//     A.zeros();
//     B.zeros();
//     for (arma::uword i = 0; i < N; i++) {
//       for (arma::uword t = 0; t < Ti(i); t++) {
//         log_py.col(t) = log_B.col(obs(t, i));
//       }
//       univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
//       univariate_backward_nhmm(log_beta, log_A, log_py.cols(0, Ti(i) - 1));
//       double ll_i = logSumExp(log_alpha.col(Ti(i) - 1));
//       ll_new += ll_i;
//       
//       // update parameters even if already converged
//       // Pi
//       Pi += arma::exp(log_alpha.col(0) + log_beta.col(0) - ll_i);
//       // A
//       for (arma::uword j = 0; j < S; j++) {
//         for (arma::uword k = 0; k < S; k++) {
//           for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
//             tmp_T(t) = log_alpha(k, t) + log_A(k, j) + log_beta(j, t + 1) + log_py(j, t + 1);
//           }
//           A(k, j) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 2) - ll_i));
//         }
//       }
//       // B
//       for (arma::uword m = 0; m < M; m++) {
//         for (arma::uword s = 0; s < S; s++) {
//           tmp_T.fill(-arma::datum::inf);
//           for (arma::uword t = 0; t < Ti(i); t++) {
//             if (m == obs(t, i)) {
//               tmp_T(t) = log_alpha(s, t) + log_beta(s, t);
//             }
//           }
//           B(s, m) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 1) - ll_i));
//         }
//       }
//     }
//     log_Pi = arma::log(Pi / arma::accu(Pi));
//     A.each_col() /= sum(A, 1);
//     log_A = arma::log(A);
//     B.cols(0, M - 1).each_col() /= sum(B.cols(0, M - 1), 1);
//     log_B.cols(0, M - 1) = arma::log(B.cols(0, M - 1));
//     
//     relative_change = (ll_new - ll) / (std::abs(ll) + 1e-8);
//     absolute_change = (ll_new - ll) / n_obs;
//     ll = ll_new;
//   }
//   Pi = arma::exp(log_Pi);
//   A = arma::exp(log_A);
//   B = arma::exp(log_B);
//   // Final log-likelihood
//   ll_new = 0;
//   for (arma::uword i = 0; i < N; i++) {
//     for (arma::uword t = 0; t < Ti(i); t++) {
//       log_py.col(t) = log_B.col(obs(t, i));
//     }
//     univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
//     ll_new += logSumExp(log_alpha.col(Ti(i) - 1));
//   }
//   return Rcpp::List::create(
//     Rcpp::Named("initial_probs") = Rcpp::wrap(Pi),
//     Rcpp::Named("transition_probs") = Rcpp::wrap(A),
//     Rcpp::Named("emission_probs") = Rcpp::wrap(B.cols(0, M - 1)),
//     Rcpp::Named("logLik") = ll_new,
//     Rcpp::Named("iterations") = iter,
//     Rcpp::Named("relative_change") = relative_change,
//     Rcpp::Named("absolute_change") = absolute_change
//   );
// }
// // // // [[Rcpp::export]]
// // // Rcpp::List EM_nhmm_multichannel(
// // //     const arma::mat& eta_pi, const arma::cube& eta_A, 
// // //     const arma::field<arma::cube>& eta_B, const arma::ucube& obs, 
// // //     const arma::uvec& M, const arma::uvec& Ti, arma::uword itermax, 
// // //     double tol) {
// // //   
// // //   arma::uword T = obs.n_rows;
// // //   arma::uword N = obs.n_cols;
// // //   arma::uword S = eta_A.n_slices;
// // //   arma::uword C = M.n_elem;
// // //   arma::uword maxM = arma::max(M);
// // // EM-algorithm begins
// // // 
// // // double change = tol + 1.0;
// // // arma::uword iter = 0;
// // // double ll_new = 0;
// // // double ll = -1e150;
// // // // initial values to probabilities
// // // arma::vec log_Pi = arma::log(softmax(eta_to_gamma(eta_pi)));
// // // arma::mat log_A(S, S);
// // // arma::mat log_B(S, M + 1);
// // // log_B.col(M).zeros();
// // // arma::cube gamma_A = eta_to_gamma(eta_A);
// // // arma::cube gamma_B = eta_to_gamma(eta_B);
// // // for (arma::uword s = 0; s < S; s++) {
// // //   log_A.row(s) = arma::log((softmax(gamma_A.slice(s)).t()));
// // //   log_B.row(s).cols(0, M - 1) = arma::log(softmax(gamma_B.slice(s)).t());
// // // }
// // // arma::vec Pi(S);
// // // arma::mat A(S, S);
// // // arma::mat B(S, M);
// // // arma::mat log_alpha(S, T);
// // // arma::mat log_beta(S, T);
// // // arma::mat log_py(S, T);
// // // arma::vec tmp_T(T);
// // // // ll = 0;
// // // // for (arma::uword i = 0; i < N; i++) {
// // // //   for (arma::uword t = 0; t < Ti(i); t++) {
// // // //     log_py.col(t) = log_B.col(obs(t, i));
// // // //   }
// // // //   univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// // // //   ll += logSumExp(log_alpha.col(Ti(i) - 1));
// // // // }
// // // while ((change > tol) && (iter < itermax)) {
// // //   iter++;
// // //   ll_new = 0;
// // //   Pi.zeros();
// // //   A.zeros();
// // //   B.zeros();
// // //   for (arma::uword i = 0; i < N; i++) {
// // //     for (arma::uword t = 0; t < Ti(i); t++) {
// // //       log_py.col(t) = log_B.col(obs(t, i));
// // //     }
// // //     univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// // //     univariate_backward_nhmm(log_beta, log_A, log_py.cols(0, Ti(i) - 1));
// // //     double ll_i = logSumExp(log_alpha.col(Ti(i) - 1));
// // //     ll_new += ll_i;
// // //     
// // //     // Pi
// // //     Pi += arma::exp(log_alpha.col(0) + log_beta.col(0) - ll_i);
// // //     // A
// // //     for (arma::uword j = 0; j < S; j++) {
// // //       for (arma::uword k = 0; k < S; k++) {
// // //         for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
// // //           tmp_T(t) = log_alpha(k, t) + log_A(k, j) + log_beta(j, t + 1) + log_py(j, t + 1);
// // //         }
// // //         A(k, j) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 2) - ll_i));
// // //       }
// // //     }
// // //     // B
// // //     for (arma::uword m = 0; m < M; m++) {
// // //       for (arma::uword s = 0; s < S; s++) {
// // //         tmp_T.fill(-arma::datum::inf);
// // //         for (arma::uword t = 0; t < Ti(i); t++) {
// // //           if (m == obs(t, i)) {
// // //             tmp_T(t) = log_alpha(s, t) + log_beta(s, t);
// // //           }
// // //         }
// // //         B(s, m) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 1) - ll_i));
// // //       }
// // //     }
// // //   }
// // //   log_Pi = arma::log(Pi / arma::accu(Pi));
// // //   A.each_col() /= sum(A, 1);
// // //   log_A = arma::log(A);
// // //   B.cols(0, M - 1).each_col() /= sum(B.cols(0, M - 1), 1);
// // //   log_B.cols(0, M - 1) = arma::log(B.cols(0, M - 1));
// // //   // ll_new = 0;
// // //   // for (arma::uword i = 0; i < N; i++) {
// // //   //   for (arma::uword t = 0; t < Ti(i); t++) {
// // //   //     log_py.col(t) = log_B.col(obs(t, i));
// // //   //   }
// // //   //   univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// // //   //   ll_new += logSumExp(log_alpha.col(Ti(i) - 1));
// // //   // }
// // //   
// // //   change = (ll_new - ll) / (std::abs(ll) + 0.1);
// // //   ll = ll_new;
// // // }
// // // Pi = arma::exp(log_Pi);
// // // A = arma::exp(log_A);
// // // B = arma::exp(log_B);
// // // // should compute the final log-likelihood using these values, 
// // // // but not interested in that here
// // // return Rcpp::List::create(
// // //   Rcpp::Named("initial_probs") = Rcpp::wrap(Pi),
// // //   Rcpp::Named("transition_probs") = Rcpp::wrap(A), 
// // //   Rcpp::Named("emission_probs") = Rcpp::wrap(B.cols(0, M - 1)),
// // //   Rcpp::Named("logLik") = ll, 
// // //   Rcpp::Named("iterations") = iter, 
// // //   Rcpp::Named("change") = change
// // // );
// // // }
