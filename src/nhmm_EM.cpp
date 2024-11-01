// // EM algorithm for NHMMs
// 
// #include "softmax.h"
// #include "forward_nhmm.h"
// #include "backward_nhmm.h"
// #include "get_parameters.h"
// #include "nhmm_mstep.h"
// #include "nhmm_sc.h"
// #include "nhmm_mc.h"
// #include "mnhmm_sc.h"
// #include "mnhmm_mc.h"

// // [[Rcpp::export]]
// Rcpp::List EM_LBFGS_nhmm_singlechannel(
//     const arma::mat& eta_pi, const arma::mat& X_pi,
//     const arma::cube& eta_A, const arma::cube& X_A,
//     const arma::cube& eta_B, const arma::cube& X_B,
//     const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B, 
//     const bool tv_A, const bool tv_B, const arma::uvec& Ti, double n_obs, 
//     arma::uword maxeval, double ftol_abs, double ftol_rel, double xtol_abs, 
//     double xtol_rel) {
//   
//   nhmm_sc model(
// eta_A.n_slices, X_pi, X_A, X_B, Ti,
// iv_pi, iv_A, iv_B, tv_A, tv_B, obs, eta_pi, eta_A, eta_B
// );

//   arma::mat log_alpha(model.S, model.T);
//   arma::mat log_beta(model.S, model.T);
//   // EM-algorithm begins
//   
//   double relative_change = ftol_rel + 1.0;
//   double absolute_change = ftol_abs + 1.0;
//   arma::uword iter = 0;
//   double ll_new = 0;
//   double ll = -1e150;
//   arma::vec tmp_T(T);
//   arma::vec E_Pi(S);
//   arma::cube E_A(S, S, T);
//   arma::cube E_B(S, M, T);
//   while ((relative_change > ftol_rel) && (absolute_change > ftol_abs) && (iter < maxeval)) {
//     iter++;
//     ll_new = 0;
//     E_Pi.zeros();
//     E_A.zeros();
//     E_B.zeros();
//     for (arma::uword i = 0; i < N; i++) {
//       
// if (model.iv_pi || i == 0) {
//   model.update_pi(i);
// }
// if (model.iv_A || i == 0) {
//   model.update_A(i);
// }
// if (model.iv_B || i == 0) {
//   model.update_B(i);
// }
// model.update_log_py(i);
// univariate_forward_nhmm(
//   log_alpha, model.log_Pi, model.log_A, 
//   model.log_py.cols(0, model.Ti(i) - 1)
// );
// univariate_backward_nhmm(
//   log_beta, model.log_A, model.log_py.cols(0, model.Ti(i) - 1)
// );
//       
//       double ll_i = logSumExp(log_alpha.col(model.Ti(i) - 1));
//       ll_new += ll_i;
//       
//       // update parameters once more even if already converged
//       // Pi
//       E_Pi += arma::exp(log_alpha.col(0) + log_beta.col(0) - ll_i);
//       // A
//       for (arma::uword j = 0; j < model.S; j++) {
//         for (arma::uword k = 0; k < model.S; k++) {
//           for (arma::uword t = 0; t < (model.Ti(i) - 1); t++) {
//             E_A(k, j, t) += log_alpha(k, t) + model.log_A(k, j, t) + log_beta(j, t + 1) + model.log_py(j, t + 1) - ll_i;
//           }
//         }
//       }
//       // B
//       for (arma::uword m = 0; m < model.M; m++) {
//         for (arma::uword s = 0; s < model.S; s++) {
//           tmp_T.fill(-arma::datum::inf);
//           for (arma::uword t = 0; t < model.Ti(i); t++) {
//             if (m == model.obs(t, i)) {
//              E_B(s, m, t) += log_alpha(s, t) + log_beta(s, t) - ll_i;
//             }
//           }
//         }
//       }
//     }
//     // Minimize obj(E_pi, E_A, E_B, eta_pi, eta_A, eta_B, X_pi, X_A, X_B)
//     // with respect to eta_pi, eta_A, eta_B
//     nhmm_mstep(model.eta_pi, model.eta_A, model.eta_B, model.X_pi, model.X_A, model.X_B, E_Pi, E_A, E_B);
//     relative_change = (ll_new - ll) / (std::abs(ll) + 1e-8);
//     absolute_change = (ll_new - ll) / n_obs;
//     ll = ll_new;
//   }
//   // Final log-likelihood
//   ll_new = 0;
//   for (arma::uword i = 0; i < N; i++) {
//     model.update_log_py(i);
//     univariate_forward_nhmm(log_alpha, model.log_Pi, model.log_A, model.log_py.cols(0, model.Ti(i) - 1));
//     ll_new += logSumExp(log_alpha.col(Ti(i) - 1));
//   }
//   return Rcpp::List::create(
//     Rcpp::Named("initial_probs") = Rcpp::wrap(model.Pi),
//     Rcpp::Named("transition_probs") = Rcpp::wrap(model.A), 
//     Rcpp::Named("emission_probs") = Rcpp::wrap(model.B.cols(0, model.M - 1)),
//     Rcpp::Named("logLik") = ll_new, 
//     Rcpp::Named("iterations") = iter, 
//     Rcpp::Named("relative_change") = relative_change, 
//     Rcpp::Named("absolute_change") = absolute_change
//   );
// }
// // // [[Rcpp::export]]
// // Rcpp::List EM_nhmm_multichannel(
// //     const arma::mat& eta_pi, const arma::cube& eta_A, 
// //     const arma::field<arma::cube>& eta_B, const arma::ucube& obs, 
// //     const arma::uvec& M, const arma::uvec& Ti, arma::uword itermax, 
// //     double tol) {
// //   
// //   arma::uword T = obs.n_rows;
// //   arma::uword N = obs.n_cols;
// //   arma::uword S = eta_A.n_slices;
// //   arma::uword C = M.n_elem;
// //   arma::uword maxM = arma::max(M);
// // EM-algorithm begins
// // 
// // double change = tol + 1.0;
// // arma::uword iter = 0;
// // double ll_new = 0;
// // double ll = -1e150;
// // // initial values to probabilities
// // arma::vec log_Pi = arma::log(softmax(eta_to_gamma(eta_pi)));
// // arma::mat log_A(S, S);
// // arma::mat log_B(S, M + 1);
// // log_B.col(M).zeros();
// // arma::cube gamma_A = eta_to_gamma(eta_A);
// // arma::cube gamma_B = eta_to_gamma(eta_B);
// // for (arma::uword s = 0; s < S; s++) {
// //   log_A.row(s) = arma::log((softmax(gamma_A.slice(s)).t()));
// //   log_B.row(s).cols(0, M - 1) = arma::log(softmax(gamma_B.slice(s)).t());
// // }
// // arma::vec Pi(S);
// // arma::mat A(S, S);
// // arma::mat B(S, M);
// // arma::mat log_alpha(S, T);
// // arma::mat log_beta(S, T);
// // arma::mat log_py(S, T);
// // arma::vec tmp_T(T);
// // // ll = 0;
// // // for (arma::uword i = 0; i < N; i++) {
// // //   for (arma::uword t = 0; t < Ti(i); t++) {
// // //     log_py.col(t) = log_B.col(obs(t, i));
// // //   }
// // //   univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// // //   ll += logSumExp(log_alpha.col(Ti(i) - 1));
// // // }
// // while ((change > tol) && (iter < itermax)) {
// //   iter++;
// //   ll_new = 0;
// //   Pi.zeros();
// //   A.zeros();
// //   B.zeros();
// //   for (arma::uword i = 0; i < N; i++) {
// //     for (arma::uword t = 0; t < Ti(i); t++) {
// //       log_py.col(t) = log_B.col(obs(t, i));
// //     }
// //     univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// //     univariate_backward_nhmm(log_beta, log_A, log_py.cols(0, Ti(i) - 1));
// //     double ll_i = logSumExp(log_alpha.col(Ti(i) - 1));
// //     ll_new += ll_i;
// //     
// //     // Pi
// //     Pi += arma::exp(log_alpha.col(0) + log_beta.col(0) - ll_i);
// //     // A
// //     for (arma::uword j = 0; j < S; j++) {
// //       for (arma::uword k = 0; k < S; k++) {
// //         for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
// //           tmp_T(t) = log_alpha(k, t) + log_A(k, j) + log_beta(j, t + 1) + log_py(j, t + 1);
// //         }
// //         A(k, j) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 2) - ll_i));
// //       }
// //     }
// //     // B
// //     for (arma::uword m = 0; m < M; m++) {
// //       for (arma::uword s = 0; s < S; s++) {
// //         tmp_T.fill(-arma::datum::inf);
// //         for (arma::uword t = 0; t < Ti(i); t++) {
// //           if (m == obs(t, i)) {
// //             tmp_T(t) = log_alpha(s, t) + log_beta(s, t);
// //           }
// //         }
// //         B(s, m) += exp(logSumExp(tmp_T.rows(0, Ti(i) - 1) - ll_i));
// //       }
// //     }
// //   }
// //   log_Pi = arma::log(Pi / arma::accu(Pi));
// //   A.each_col() /= sum(A, 1);
// //   log_A = arma::log(A);
// //   B.cols(0, M - 1).each_col() /= sum(B.cols(0, M - 1), 1);
// //   log_B.cols(0, M - 1) = arma::log(B.cols(0, M - 1));
// //   // ll_new = 0;
// //   // for (arma::uword i = 0; i < N; i++) {
// //   //   for (arma::uword t = 0; t < Ti(i); t++) {
// //   //     log_py.col(t) = log_B.col(obs(t, i));
// //   //   }
// //   //   univariate_forward_nhmm(log_alpha, log_Pi, log_A, log_py.cols(0, Ti(i) - 1));
// //   //   ll_new += logSumExp(log_alpha.col(Ti(i) - 1));
// //   // }
// //   
// //   change = (ll_new - ll) / (std::abs(ll) + 0.1);
// //   ll = ll_new;
// // }
// // Pi = arma::exp(log_Pi);
// // A = arma::exp(log_A);
// // B = arma::exp(log_B);
// // // should compute the final log-likelihood using these values, 
// // // but not interested in that here
// // return Rcpp::List::create(
// //   Rcpp::Named("initial_probs") = Rcpp::wrap(Pi),
// //   Rcpp::Named("transition_probs") = Rcpp::wrap(A), 
// //   Rcpp::Named("emission_probs") = Rcpp::wrap(B.cols(0, M - 1)),
// //   Rcpp::Named("logLik") = ll, 
// //   Rcpp::Named("iterations") = iter, 
// //   Rcpp::Named("change") = change
// // );
// // }
