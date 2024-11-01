// #include "eta_to_gamma.h"
// #include "get_parameters.h"
// #include "logsumexp.h"
// #include "nhmm.h"
// #include "sum_to_zero.h"
// #include <nloptrAPI.h>
// 
// struct nhmm_opt_data {
//   nhmm model;
//   arma::vec E_Pi;
//   arma::mat E_A;
//   arma::mat E_B;
//   nhmm_opt_data(const nhmm& model_, const arma::vec& E_Pi_, 
//                 const arma::vec& E_A_, const arma::vec& E_B_): model(model_), 
//                 E_Pi(E_Pi_), E_A(E_A_), E_B(E_B_) {}
// };
// // Global counter for function evaluations
// static int iter = 0;
// 
// // Define actual objective function
// double objective_pi(const arma::vec& x, arma::vec& grad, void *data) {
//   iter++;
//   // cast generic void* pointer data to a pointer to a nhmm_opt_data type
//   nhmm_opt_data* opt_data = static_cast<nhmm_opt_data*>(data);
//   // model is a pointer to a model of opt_data
//   nhmm* model = &(opt_data->model);
//   model->eta_pi = arma::mat(x.memptr(), model->S - 1, model->K_pi);
//   model->gamma_pi = sum_to_zero(model->eta_pi, model->Qs);
//   arma::vec Pi(model->S);
//   double value = 0;
//   arma::mat Qt = (model->Qs).t();
//   for (arma::uword i = 0; i < model->N; i++) {
//     if (model->iv_pi || i == 0) {
//       Pi = get_pi(model->gamma_pi, model->X_pi.col(i));
//     }
//     value -= arma::as_scalar(opt_data->E_Pi.t() * model->gamma_pi * model->X_pi.col(i) - sum(opt_data->E_Pi) * logSumExp(model->gamma_pi * model->X_pi.col(i)));
//     grad += Qt * (opt_data->E_Pi - sum(opt_data->E_Pi) * Pi) * model->X_pi.col(i).t();
//   }
//   return value;
// }
// 
// // Wrapper function for NLopt to interface with Armadillo
// double c_objective_pi(unsigned n, const double *x, double *grad, void *data) {
//   // Wrap raw pointer `x` in an Armadillo vector (read/write access)
//   arma::vec x_vec(const_cast<double *>(x), n, false, true); 
//   // Wrap the gradient pointer `grad` in an Armadillo vector
//   arma::vec grad_vec(grad, n, false, true);  
//   // Call the Armadillo-based objective function
//   return objective_pi(x_vec, grad_vec, data);  
// }
// 
// // // Wrapper function for NLopt to interface with Armadillo
// // double c_objective_A(unsigned n, const double *x, double *grad, void *data) {
// //   // Wrap raw pointer `x` in an Armadillo vector (read/write access)
// //   arma::vec x_vec(const_cast<double *>(x), n, false, true); 
// //   
// //   // Wrap the gradient pointer `grad` in an Armadillo vector
// //   arma::vec grad_vec(grad, n, false, true);  
// //   
// //   // Call the Armadillo-based function
// //   return objective_A(x_vec, grad_vec, data);  
// // }
// // 
// // // Wrapper function for NLopt to interface with Armadillo
// // double c_objective_B(unsigned n, const double *x, double *grad, void *data) {
// //   // Wrap raw pointer `x` in an Armadillo vector (read/write access)
// //   arma::vec x_vec(const_cast<double *>(x), n, false, true); 
// //   
// //   // Wrap the gradient pointer `grad` in an Armadillo vector
// //   arma::vec grad_vec(grad, n, false, true);  
// //   
// //   // Call the Armadillo-based function
// //   return objective_B(x_vec, grad_vec, data);  
// // }
// 
// 
// // double objective_A(const arma::vec& x, arma::vec& grad, void *data) {
// //   iter++;
// //   nhmm* model = static_cast<nhmm*>(data.model);
// //   arma::vec E_A = static_cast<arma::vec*>(data.E_A);
// //   unsigned state = static_cast<arma::uword>(data.state);
// //   arma::mat eta_A = arma::mat(x.memptr(), model->S - 1, model->K_A);
// //   arma::mat gamma_A = sum_to_zero(model->eta_A, model->Qs);
// //   arma::vec A(model->S);
// //   double value = 0;
// //   arma::mat Qt = (model->Qs).t();
// //   for (arma::uword i = 0; i < model->N; i++) {
// //     if (model->iv_A || i == 0) {
// //       A = get_pi(model->gamma_A, model->X_A.col(i));
// //     }
// //     value -= (E_Pi.t() * gamma_pi * X_pi.col(i) - sum(E_Pi) * logSumExp(model->gamma_pi * model->X_pi.col(i)));
// //     grad += Qt * (E_Pi - sum(E_Pi) * Pi) * model->X.col(i).t();
// //   }
// //   return value;
// // }
// // double objective_pi(const arma::vec& x, arma::vec& grad, void *data) {
// //   iter++;
// //   nhmm* model = static_cast<nhmm*>(data);
// //   model->eta_pi = arma::mat(x.memptr(), model->S - 1, model->K_pi);
// //   model->gamma_pi = sum_to_zero(model->eta_pi, model->Qs);
// //   arma::vec Pi(model->S);
// //   double value = 0;
// //   arma::mat Qt = (model->Qs).t();
// //   for (arma::uword i = 0; i < model->N; i++) {
// //     if (model->iv_pi || i == 0) {
// //       Pi = get_pi(model->gamma_pi, model->X_pi.col(i));
// //     }
// //     value -= (E_Pi.t() * gamma_pi * X_pi.col(i) - sum(E_Pi) * logSumExp(model->gamma_pi * model->X_pi.col(i)));
// //     grad += Qt * (E_Pi - sum(E_Pi) * Pi) * model->X.col(i).t();
// //   }
// //   return value;
// // }
// void nhmm::mstep(
//     const arma::vec E_Pi, const arma::cube E_A, const arma::cube E_B,
//     const double xtol_abs, const double ftol_abs, const double xtol_rel, 
//     const double ftol_rel, arma::uword maxeval) {
//   
//   nhmm_opt_data opt_data(model, E_Pi, E_A, E_B); 
//   // pi
//   nlopt_opt opt_pi = nlopt_create(NLOPT_LD_LBFGS, model.np_pi);
//   nlopt_set_min_objective(opt_pi, c_objective_pi, &opt_data);
//   nlopt_set_xtol_abs1(opt_pi, xtol_abs);
//   nlopt_set_ftol_abs(opt_pi, ftol_abs);
//   nlopt_set_xtol_rel(opt_pi, xtol_rel);
//   nlopt_set_ftol_rel(opt_pi, ftol_rel);
//   nlopt_set_maxeval(opt_pi, maxeval);
//   arma::vec x = arma::vectorise(model.eta_pi); 
//   double minf;
//   iter = 0; // Reset counter
//   int status = nlopt_optimize(opt_pi, x.memptr(), &minf);
//   nlopt_destroy(opt_pi);
//   
//   // // A
//   // nlopt_opt opt_A = nlopt_create(NLOPT_LD_LBFGS, model.np_pi);
//   // nlopt_set_min_objective(opt_A, c_objective_A, &model);
//   // nlopt_set_xtol_abs(opt_A, xtol_abs);
//   // nlopt_set_ftol_abs(opt_A, ftol_abs);
//   // nlopt_set_xtol_rel(opt_A, xtol_rel);
//   // nlopt_set_ftol_rel(opt_A, ftol_rel);
//   // nlopt_set_maxeval(opt_A, maxeval);
//   // arma::vec x = arma::vectorise(model.eta_A); 
//   // double minf;
//   // iter = 0; // Reset counter
//   // int status = nlopt_optimize(opt_A, x.memptr(), &minf);
//   // nlopt_destroy(opt_A);
//   // 
//   // // B
//   // nlopt_opt opt_B = nlopt_create(NLOPT_LD_LBFGS, model.np_pi);
//   // nlopt_set_min_objective(opt_B, c_objective_B, &model);
//   // nlopt_set_xtol_abs(opt_B, xtol_abs);
//   // nlopt_set_ftol_abs(opt_B, ftol_abs);
//   // nlopt_set_xtol_rel(opt_B, xtol_rel);
//   // nlopt_set_ftol_rel(opt_B, ftol_rel);
//   // nlopt_set_maxeval(opt_B, maxeval);
//   // arma::vec x = arma::vectorise(model.eta_B); 
//   // double minf;
//   // iter = 0; // Reset counter
//   // int status = nlopt_optimize(opt_B, x.memptr(), &minf);
//   // nlopt_destroy(opt_B);
//   // 
// }
// // 
// // double nhmm_sc_obj(
// //     arma::vec& grad,
// //     const arma::mat& Qs, const arma::mat& Qm,
// //     arma::mat& eta_pi, const arma::mat& X_pi,
// //     arma::cube& eta_A, const arma::cube& X_A,
// //     arma::cube& eta_B, const arma::cube& X_B,
// //     const arma::umat& obs, const bool iv_pi, const bool iv_A, const bool iv_B,
// //     const bool tv_A, const bool tv_B, const arma::uvec& Ti,
// //     const arma::vec E_Pi, const arma::cube E_A, const arma::cube E_B) {
// //   
// //   arma::uword N = X_A.n_slices;
// //   arma::uword T = X_A.n_cols;
// //   arma::uword S = eta_A.n_slices;
// //   arma::uword M = eta_B.n_rows + 1;
// //   
// //   arma::mat gamma_pi = eta_to_gamma(eta_pi, Qs);
// //   arma::cube gamma_A = eta_to_gamma(eta_A, Qs);
// //   arma::cube gamma_B = eta_to_gamma(eta_B, Qm);
// //   arma::vec log_Pi(S);
// //   arma::cube log_A(S, S, T);
// //   arma::cube log_B(S, M, T);
// //   arma::vec Pi(S);
// //   arma::cube A(S, S, T);
// //   arma::cube B(S, M, T);
// //   double nll = 0;
// //   for (arma::uword i = 0; i < N; i++) {
// //     
// //     if (iv_pi || i == 0) {
// //       Pi = get_log_pi(gamma_pi, X_pi.col(i));
// //       log_Pi = arma::log(Pi);
// //     }
// //     if (iv_A || i == 0) {
// //       A = get_A(gamma_A, X_A.slice(i), tv_A);
// //       log_A = arma::log(A);
// //     }
// //     if (iv_B || i == 0) {
// //       B = get_B(gamma_B, X_B.slice(i), tv_B, false);
// //       log_B = arma::log(B);
// //     }
// //     for (arma::uword s = 0; s < S; ++s) {
// //       double pi_s = arma::dot(gamma_pi.col(s), X_pi.row(s));
// //       nll -= E_Pi(s) * std::log(pi_s);
// //       for (arma::uword j = 0; j < gamma_pi.n_rows; ++j) {
// //         grad_pi(j, s) -= E_Pi(s) * (X_pi(j, s) - pi_s * gamma_pi(j, s)) / pi_s;
// //       }
// //     }
// //     for (arma::uword s = 0; s < S; ++s) {
// //       for (arma::uword s_next = 0; s_next < S; ++s_next) {
// //         for (arma::uword t = 0; t < E_A.n_slices; ++t) {
// //           double A_ss = arma::dot(gamma_A.slice(s).col(s_next), X_A.slice(t).row(s));
// //           nll -= E_A(s, s_next, t) * std::log(A_ss);
// //           for (arma::uword j = 0; j < gamma_A.n_rows; ++j) {
// //             grad_A(j, s, s_next) -= E_A(s, s_next, t) * (X_A(j, t) - A_ss * gamma_A(j, s)) / A_ss;
// //           }
// //         }
// //       }
// //     }
// //     
// //     // Emission Matrix Contribution to NLL and Gradients
// //     for (arma::uword s = 0; s < S; ++s) {
// //       for (arma::uword m = 0; m < M; ++m) {
// //         for (arma::uword t = 0; t < E_B.n_slices; ++t) {
// //           double B_sm = arma::dot(gamma_B.slice(s).col(m), X_B.slice(t).row(s));
// //           nll -= E_B(s, m, t) * std::log(B_sm);
// //           for (arma::uword j = 0; j < gamma_B.n_rows; ++j) {
// //             grad_B(j, s, m) -= E_B(s, m, t) * (X_B(j, t) - B_sm * gamma_B(j, s)) / B_sm;
// //           }
// //         }
// //       }
// //     }
// //     
// //     arma::uword N = X_A.n_slices;
// //     arma::uword T = X_A.n_cols;
// //     arma::uword S = eta_A.n_slices;
// //     arma::uword M = eta_B.n_rows + 1;
// //     arma::vec loglik(N);
// //     arma::mat log_alpha(S, T);
// //     arma::mat log_beta(S, T);
// //     arma::mat log_py(S, T);
// //     
// //     arma::vec Pi(S);
// //     arma::cube A(S, S, T);
// //     arma::cube B(S, M + 1, T);
// //     arma::vec log_Pi(S);
// //     arma::cube log_A(S, S, T);
// //     arma::cube log_B(S, M + 1, T);
// //     
// //     arma::mat grad_pi(S - 1, X_pi.n_rows, arma::fill::zeros);
// //     arma::cube grad_A(S - 1, X_A.n_rows, S, arma::fill::zeros);
// //     arma::cube grad_B(M - 1, X_B.n_rows, S, arma::fill::zeros);
// //     
// //     arma::mat gamma_pi = eta_to_gamma(eta_pi, Qs.t());
// //     arma::cube gamma_A = eta_to_gamma(eta_A, Qs.t());
// //     arma::cube gamma_B = eta_to_gamma(eta_B, Qm.t());
// //     for (arma::uword i = 0; i < N; i++) {
// //       if (iv_pi || i == 0) {
// //         Pi = get_pi(gamma_pi, X_pi.col(i));
// //         log_Pi = arma::log(Pi);
// //       }
// //       if (iv_A || i == 0) {
// //         A = get_A(gamma_A, X_A.slice(i), tv_A);
// //         log_A = arma::log(A);
// //       }
// //       if (iv_B || i == 0) {
// //         B = get_B(gamma_B, X_B.slice(i), true, tv_B);
// //         log_B = arma::log(B);
// //       }
// //       for (arma::uword t = 0; t < Ti(i); t++) {
// //         log_py.col(t) = log_B.slice(t).col(obs(t, i));
// //       }
// //       log_alpha = univariate_forward_nhmm(log_Pi, log_A, log_py);
// //       log_beta = univariate_backward_nhmm(log_A, log_py);
// //       double ll = logSumExp(log_alpha.col(T - 1));
// //       if (!std::isfinite(ll)) {
// //         grad.fill(arma::datum::inf);
// //         return arma::datum::inf;
// //       }
// //       loglik(i) = ll;
// //       // gradient wrt gamma_pi
// //       grad_pi += gradient_wrt_pi(Qs, log_py, log_beta, ll, Pi, X_pi, i);
// //       // gradient wrt gamma_A
// //       for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
// //         for (arma::uword s = 0; s < S; s++) {
// //           grad_A.slice(s) += gradient_wrt_A(Qs, log_py, log_alpha, log_beta, ll, A, X_A, i, t, s);
// //         }
// //       }
// //       // gradient wrt gamma_B
// //       for (arma::uword s = 0; s < S; s++) {
// //         if (obs(0, i) < M) {
// //           grad_B.slice(s) += gradient_wrt_B_t0(Qm, obs, log_Pi, log_beta, ll, B, X_B, i, s);
// //         }
// //         for (arma::uword t = 0; t < (Ti(i) - 1); t++) {
// //           if (obs(t + 1, i) < M) {
// //             grad_B.slice(s) += gradient_wrt_B(Qm, obs, log_alpha, log_beta, ll, log_A, B, X_B, i, s, t);
// //           }
// //         }
// //       }
// //     }
// //     size_t n_pi = (S - 1) * X_pi.n_rows;
// //     size_t n_A = (S - 1) * X_A.n_rows * S;
// //     size_t n_B = (M - 1) * X_B.n_rows * S;
// //     size_t idx = 0;
// //     grad.subvec(idx, n_pi - 1) = -arma::vectorise(grad_pi);
// //     idx += n_pi;
// //     grad.subvec(idx, idx + n_A - 1) = -arma::vectorise(grad_A);
// //     idx += n_A;
// //     grad.subvec(idx, idx + n_B - 1) = -arma::vectorise(grad_B);
// //     return -sum(loglik);
// //   }
