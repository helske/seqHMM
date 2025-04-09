#ifndef NHMM_H
#define NHMM_H

#include "config.h"
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"
#include <nloptrAPI.h>

struct nhmm {
  const arma::ucube& obs;
  const arma::uvec& Ti;
  const arma::uvec& M;
  const arma::uword N;
  const arma::uword T;
  const arma::uword C;
  const arma::uword S;
  const arma::mat& X_pi;
  const arma::cube& X_A;
  const arma::field<arma::cube>& X_B;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const arma::uvec& icpt_only_B;
  const bool iv_A;
  const arma::uvec& iv_B;
  const bool tv_A; 
  const arma::uvec& tv_B;
  const arma::mat Qs;
  const arma::field<arma::mat> Qm;
  //coefficients //
  arma::mat eta_pi;
  arma::cube eta_A;
  arma::field<arma::cube> eta_B;  
  arma::mat gamma_pi;
  arma::cube gamma_A;
  arma::field<arma::cube> gamma_B;
  
  // pi, A, B, and log_p(y) of _one_ id we are currently working with
  arma::mat log_py;
  arma::vec pi;
  arma::vec log_pi;
  arma::cube A;
  arma::cube log_A;
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  
  // excepted counts for EM algorithm
  arma::mat E_pi;
  arma::field<arma::cube> E_A;
  arma::field<arma::cube> E_B;
  arma::uword current_s;
  arma::uword current_c; 
  
  const arma::uword n_obs;
  const double lambda; // regularization
  const double maxval;
  const double minval;
  // for EM
  int mstep_iter = 0;
  int mstep_return_code = 0;
  nlopt_opt opt_pi = nullptr;
  nlopt_opt opt_A = nullptr;
  std::vector<nlopt_opt> opt_B;
  
  nhmm(
    const arma::ucube& obs_,
    const arma::uvec& Ti_,
    const arma::uvec& M_,
    const arma::mat& X_pi_,
    const arma::cube& X_A_,
    const arma::field<arma::cube>& X_B_,
    const bool icpt_only_pi_,
    const bool icpt_only_A_,
    const arma::uvec& icpt_only_B_,
    const bool iv_A_,
    const arma::uvec& iv_B_,
    const bool tv_A_,
    const arma::uvec& tv_B_,
    const arma::mat& eta_pi_,
    const arma::cube& eta_A_,
    const arma::field<arma::cube>& eta_B_,
    const double lambda_ = 0,
    double maxval_ = arma::datum::inf,
    double minval_ = -1.0)
    :
    obs(obs_),
    Ti(Ti_),
    M(M_),
    N(obs.n_slices),
    T(obs.n_cols),
    C(obs.n_rows),
    S(eta_A_.n_slices),
    X_pi(X_pi_),
    X_A(X_A_),
    X_B(X_B_),
    icpt_only_pi(icpt_only_pi_),
    icpt_only_A(icpt_only_A_),
    icpt_only_B(icpt_only_B_),
    iv_A(iv_A_),
    iv_B(iv_B_),
    tv_A(tv_A_),
    tv_B(tv_B_),
    Qs(create_Q(S)),
    Qm(create_Q(M_)),
    eta_pi(eta_pi_),
    eta_A(eta_A_),
    eta_B(eta_B_),
    gamma_pi(eta_to_gamma(eta_pi, Qs)),
    gamma_A(eta_to_gamma(eta_A, Qs)),
    gamma_B(eta_to_gamma(eta_B, Qm)),
    log_py(S, T),
    pi(S),
    log_pi(S),
    A(S, S, T),
    log_A(S, S, T),
    B(C),
    log_B(C),
    E_pi(S, N),
    E_A(S),
    E_B(C),
    current_s(0),
    current_c(0),
    n_obs(arma::accu(Ti_)),
    lambda(lambda_),
    maxval(maxval_),
    minval(
      (minval_ < 0)  ? std::pow(arma::datum::eps, 2.0/3.0) : minval_
    ),
    opt_B(C, nullptr) {
    for (arma::uword s = 0; s < S; ++s) {
      E_A(s) = arma::cube(S, N, T, arma::fill::zeros);
    }
    for (arma::uword c = 0; c < C; ++c) {
      B(c) = arma::cube(S, M(c) + 1, T);   // B field initialization
      log_B(c) = arma::cube(S, M(c) + 1, T); // log_B field initialization
      E_B(c) = arma::cube(T, N, S);
    }
  }
  ~nhmm() {
    if (opt_pi) {
      nlopt_destroy(opt_pi);
    }
    if (opt_A) {
      nlopt_destroy(opt_A);
    }    
    for (auto& opt : opt_B) {
      if (opt) {
        nlopt_destroy(opt);  // Clean up the optimizer
        opt = nullptr;  // Prevent dangling pointer
      }
    }
  }
  void update_gamma_pi() {
    gamma_pi = eta_to_gamma(eta_pi, Qs);
  }
  void update_gamma_A() {
    gamma_A = eta_to_gamma(eta_A, Qs);
  }
  void update_gamma_B() {
    for (arma::uword c = 0; c < C; ++c) {
      gamma_B(c) = eta_to_gamma(eta_B(c), Qm(c));
    }
  }
  void update_pi(arma::uword i) {
    if (icpt_only_pi) {
      pi = softmax(gamma_pi.col(0));
    } else {
      pi = softmax(gamma_pi * X_pi.col(i));
    }
    log_pi = arma::log(pi);
  }
  void update_A(arma::uword i) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Atmp.col(s) = softmax(gamma_A.slice(s).col(0));
      }
      A.each_slice() = Atmp.t();
    } else {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); ++t) { // time
          for (arma::uword s = 0; s < S; ++s) { // from states
            Atmp.col(s) = softmax(gamma_A.slice(s) * X_A.slice(i).col(t));
          }
          A.slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A.slice(s) * X_A.slice(i).col(0));
        }
        A.each_slice() = Atmp.t();
      }
    }
    log_A = arma::log(A);
  }
  void update_B(const arma::uword i) {
    for (arma::uword c = 0; c < C; ++c) {
      if (icpt_only_B(c)) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(c).slice(s).col(0)
          );
        }
        B(c).each_slice() = Btmp.t();
      } else {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        if (tv_B(c)) {
          for (arma::uword t = 0; t < Ti(i); ++t) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(c).slice(s) * X_B(c).slice(i).col(t));
            }
            B(c).slice(t) = Btmp.t();
          }
        } else {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(c).slice(s) * X_B(c).slice(i).col(0)
            );
          }
          B(c).each_slice() = Btmp.t();
        }
      }
      log_B(c) = arma::log(B(c));
    }
  }
  void update_log_py(const arma::uword i) {
    log_py.zeros();
    for (arma::uword t = 0; t < Ti(i); ++t) {
      for (arma::uword c = 0; c < C; ++c) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
  }
  void estep_pi(const arma::uword i, const arma::vec& log_alpha, 
                const arma::vec& log_beta, const double ll) {
    E_pi.col(i) = arma::exp(log_alpha + log_beta - ll);
    // set minuscule values to zero in order to avoid numerical issues
    E_pi.col(i).clean(minval);
  }
  void estep_A(const arma::uword i, const arma::mat& log_alpha, 
               const arma::mat& log_beta, const double ll) {
    for (arma::uword k = 0; k < S; ++k) { // from
      for (arma::uword j = 0; j < S; ++j) { // to
        for (arma::uword t = 0; t < (Ti(i) - 1); ++t) { // time
          E_A(k)(j, i, t + 1) = exp(log_alpha(k, t) + log_A(k, j, t + 1) + 
            log_beta(j, t + 1) + log_py(j, t + 1) - ll);
        }
      }
      // set small values to zero in order to avoid numerical issues
      E_A(k).col(i).clean(minval);
    }
  }
  void estep_B(const arma::uword i, const arma::mat& log_alpha, 
               const arma::mat& log_beta, const double ll) {
    for (arma::uword k = 0; k < S; ++k) { // state
      for (arma::uword t = 0; t < Ti(i); ++t) { // time
        double pp = exp(log_alpha(k, t) + log_beta(k, t) - ll);
        for (arma::uword c = 0; c < C; ++c) { // channel
          if (obs(c, t, i) < M(c) && pp > minval) {
            E_B(c)(t, i, k) = pp;
          } else {
            E_B(c)(t, i, k) = 0.0;
          }
        }
      }
    }
  }
  void mstep_pi(const arma::uword print_level);
  void mstep_A(const arma::uword print_level);
  void mstep_B(const arma::uword print_level);
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
  double objective_B(const arma::vec& x, arma::vec& grad);
  static double objective_pi_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<nhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_pi(x_vec, grad_vec);
  }
  static double objective_A_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<nhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_A(x_vec, grad_vec);
  }
  static double objective_B_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<nhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  }
  Rcpp::List EM(
      const arma::uword maxeval, const double ftol_abs, const double ftol_rel, 
      const double xtol_abs, const double xtol_rel, const arma::uword print_level,
      const arma::uword maxeval_m, const double ftol_abs_m, const double ftol_rel_m, 
      const double xtol_abs_m, const double xtol_rel_m,
      const arma::uword print_level_m, const double bound
  );
  Rcpp::List mstep_error(int return_code, int iter, double relative_change, 
                         double absolute_change, double absolute_x_change, 
                         double relative_x_change) {
    if (return_code != 0) {
      return Rcpp::List::create(
        Rcpp::Named("return_code") = return_code,
        Rcpp::Named("eta_pi") = Rcpp::wrap(eta_pi),
        Rcpp::Named("eta_A") = Rcpp::wrap(eta_A),
        Rcpp::Named("eta_B") = Rcpp::wrap(eta_B),
        Rcpp::Named("logLik") = arma::datum::nan,
        Rcpp::Named("penalty_term") = arma::datum::nan,
        Rcpp::Named("iterations") = iter,
        Rcpp::Named("relative_f_change") = relative_change,
        Rcpp::Named("absolute_f_change") = absolute_change,
        Rcpp::Named("absolute_x_change") = absolute_x_change,
        Rcpp::Named("relative_x_change") = relative_x_change
      );
    }
    return Rcpp::List(); // Empty list indicates no error
  }
  void predict(arma::field<arma::cube>& obs_prob);
  void viterbi(arma::umat& q, arma::vec& logp) {
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_pi || i == 0) {
        update_pi(i);
      }
      if (iv_A || i == 0) {
        update_A(i);
      }
      if (arma::any(iv_B) || i == 0) {
        update_B(i);
      }
      update_log_py(i);
      arma::subview_col<unsigned int> subcol = q.col(i);
      logp(i) = univariate_viterbi(
        subcol,
        log_pi, 
        log_A.slices(0, Ti(i) - 1), 
        log_py.cols(0, Ti(i) - 1)
      );
    }
  }
  void forward(arma::cube& log_alpha) {
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_pi || i == 0) {
        update_pi(i);
      }
      if (iv_A || i == 0) {
        update_A(i);
      }
      if (arma::any(iv_B) || i == 0) {
        update_B(i);
      }
      update_log_py(i);
      univariate_forward(
        log_alpha.slice(i),
        log_pi,
        log_A, 
        log_py.cols(0,Ti(i) - 1)
      );
    }
  }
  
  void backward(arma::cube& log_beta) {
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_pi || i == 0) {
        update_pi(i);
      }
      if (iv_A || i == 0) {
        update_A(i);
      }
      if (arma::any(iv_B) || i == 0) {
        update_B(i);
      }
      update_log_py(i);
      univariate_backward(
        log_beta.slice(i),
        log_A, 
        log_py.cols(0, Ti(i) - 1)
      );
    }
  }
};

#endif
