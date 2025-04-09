#ifndef MNHMM_H
#define MNHMM_H

#include "config.h"
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "softmax.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"
#include "list_to_2d_field.h"
#include <nloptrAPI.h>

struct mnhmm { 
  const arma::ucube& obs;
  const arma::uvec& Ti;
  const arma::uvec& M;
  const arma::uword N;
  const arma::uword T;
  const arma::uword C;
  const arma::uword S;
  const arma::uword D;
  const arma::mat& X_pi;
  const arma::cube& X_A;
  const arma::field<arma::cube>& X_B;
  const arma::mat& X_omega;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const arma::uvec& icpt_only_B;
  const bool icpt_only_omega;
  const bool iv_A;
  const arma::uvec& iv_B;
  const bool tv_A; 
  const arma::uvec& tv_B;
  
  const arma::mat Qs;
  const arma::field<arma::mat> Qm;
  const arma::mat Qd;
  
  // coefficients //
  arma::field<arma::mat> eta_pi;
  arma::field<arma::mat> gamma_pi;
  arma::field<arma::cube> eta_A;
  arma::field<arma::cube> gamma_A;
  arma::field<arma::cube> eta_B;  
  arma::field<arma::cube> gamma_B;
  arma::mat eta_omega;
  arma::mat gamma_omega;
  
  // pi, A, B, omega, and log_p(y) of _one_ id we are currently working with
  arma::cube log_py;
  arma::field<arma::vec> pi;
  arma::field<arma::vec> log_pi;
  arma::field<arma::cube> A;
  arma::field<arma::cube> log_A;
  arma::field<arma::cube> B;
  arma::field<arma::cube> log_B;
  arma::vec omega;
  arma::vec log_omega;
  
  // excepted counts for EM algorithm
  arma::field<arma::mat> E_pi;
  arma::field<arma::cube> E_A;
  arma::field<arma::cube> E_B;
  arma::mat E_omega;
  arma::uword current_s;
  arma::uword current_c; 
  arma::uword current_d;
  
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
  nlopt_opt opt_omega = nullptr;
  mnhmm(
    const arma::ucube& obs_,
    const arma::uvec& Ti_,
    const arma::uvec& M_,
    const arma::mat& X_pi_,
    const arma::cube& X_A_,
    const arma::field<arma::cube>& X_B_,
    const arma::mat& X_omega_,
    const bool icpt_only_pi_,
    const bool icpt_only_A_,
    const arma::uvec& icpt_only_B_,
    const bool icpt_only_omega_,
    const bool iv_A_,
    const arma::uvec& iv_B_,
    const bool tv_A_,
    const arma::uvec& tv_B_,
    const arma::field<arma::mat>& eta_pi_,
    const arma::field<arma::cube>& eta_A_,
    const Rcpp::List& eta_B_,
    const arma::mat& eta_omega_,
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
    S(eta_A_(0).n_slices),
    D(eta_A_.n_elem),
    X_pi(X_pi_),
    X_A(X_A_),
    X_B(X_B_),
    X_omega(X_omega_),
    icpt_only_pi(icpt_only_pi_),
    icpt_only_A(icpt_only_A_),
    icpt_only_B(icpt_only_B_),
    icpt_only_omega(icpt_only_omega_),
    iv_A(iv_A_),
    iv_B(iv_B_),
    tv_A(tv_A_),
    tv_B(tv_B_),
    Qs(create_Q(S)),
    Qm(create_Q(M)),
    Qd(create_Q(D)),
    eta_pi(eta_pi_),
    gamma_pi(eta_to_gamma(eta_pi, Qs)),
    eta_A(eta_A_),
    gamma_A(eta_to_gamma(eta_A, Qs)),
    eta_B(list_to_2d_field(eta_B_)),
    gamma_B(eta_to_gamma(eta_B, Qm, D)),
    eta_omega(eta_omega_),
    gamma_omega(eta_to_gamma(eta_omega, Qd)),
    log_py(S, T, D),
    pi(D),
    log_pi(D),
    A(D),
    log_A(D),
    B(D, C),
    log_B(D, C),
    omega(D),
    log_omega(D),
    E_pi(D),
    E_A(D, S),
    E_B(D, C),
    E_omega(D, N),
    current_s(0),
    current_c(0),
    current_d(0),
    n_obs(arma::accu(Ti_)),
    lambda(lambda_),
    maxval(maxval_),
    minval(
      (minval_ < 0)  ? std::pow(arma::datum::eps, 2.0/3.0) : minval_
    ),
    opt_B(C, nullptr)
  {
    for (arma::uword d = 0; d < D; ++d) {
      pi(d) = arma::vec(S);
      log_pi(d) = arma::vec(S);
      A(d) = arma::cube(S, S, T);
      log_A(d) = arma::cube(S, S, T);
      E_pi(d) = arma::mat(S, N);
      for (arma::uword s = 0; s < S; ++s) {
        E_A(d, s) = arma::cube(S, N, T, arma::fill::zeros);
      }
      for (arma::uword c = 0; c < C; ++c) {
        B(d, c) = arma::cube(S, M(c) + 1, T);   // B field initialization
        log_B(d, c) = arma::cube(S, M(c) + 1, T); // log_B field initialization
        E_B(d, c) = arma::cube(T, N, S);
      }
    }
  }
  
  void update_gamma_pi() {
    for (arma::uword d = 0; d < D; ++d) {
      gamma_pi(d) = eta_to_gamma(eta_pi(d), Qs);
    }
  }
  void update_gamma_A() {
    for (arma::uword d = 0; d < D; ++d) {
      gamma_A(d) = eta_to_gamma(eta_A(d), Qs);
    }
  }
  void update_gamma_B() {
    for (arma::uword c = 0; c < C; ++c) {
      for (arma::uword d = 0; d < D; ++d) {
        gamma_B(d, c) = eta_to_gamma(eta_B(d, c), Qm(c));
      }
    }
  }
  void update_gamma_omega() {
    gamma_omega = eta_to_gamma(eta_omega, Qd);
  }
  
  void update_pi(const arma::uword i) {
    if (icpt_only_pi) {
      for (arma::uword d = 0; d < D; ++d) {
        pi(d) = softmax(gamma_pi(d).col(0));
        log_pi(d) = arma::log(pi(d));
      }
    } else {
      for (arma::uword d = 0; d < D; ++d) {
        pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
        log_pi(d) = arma::log(pi(d));
      }
    }
  }
  void update_pi(const arma::uword i, const arma::uword d) {
    if (icpt_only_pi) {
      pi(d) = softmax(gamma_pi(d).col(0));
    } else {
      pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
    }
    log_pi(d) = arma::log(pi(d));
  }
  void update_A(const arma::uword i) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword d = 0; d < D; ++d) {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s).col(0));
        }
        A(d).each_slice() = Atmp.t();
        log_A(d) = arma::log(A(d));
      }
    } else {
      for (arma::uword d = 0; d < D; ++d) {
        if (tv_A) {
          for (arma::uword t = 0; t < Ti(i); ++t) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(t));
            }
            A(d).slice(t) = Atmp.t();
          }
        } else {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(0));
          }
          A(d).each_slice() = Atmp.t();
        }
        log_A(d) = arma::log(A(d));
      }
    }
  }
  void update_A(const arma::uword i, const arma::uword d) {
    arma::mat Atmp(S, S);
    if (icpt_only_A) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Atmp.col(s) = softmax(gamma_A(d).slice(s).col(0));
      }
      A(d).each_slice() = Atmp.t();
    } else {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); ++t) { // time
          for (arma::uword s = 0; s < S; ++s) { // from states
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(t));
          }
          A(d).slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(0));
        }
        A(d).each_slice() = Atmp.t();
      }
      
    }
    log_A(d) = arma::log(A(d));
  }
  void update_B(const arma::uword i) {
    for (arma::uword c = 0; c < C; ++c) {
      if (icpt_only_B(c)) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        for (arma::uword d = 0; d < D; ++d) {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(d, c).slice(s).col(0)
            );
          }
          B(d, c).each_slice() = Btmp.t();
          log_B(d, c) = arma::log(B(d, c));
        }
      } else { 
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        if (tv_B(c)) {
          for (arma::uword d = 0; d < D; ++d) {
            for (arma::uword t = 0; t < Ti(i); ++t) { // time
              for (arma::uword s = 0; s < S; ++s) { // from states
                Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(d, c).slice(s) * X_B(c).slice(i).col(t));
              }
              B(d, c).slice(t) = Btmp.t();
            }
            log_B(d, c) = arma::log(B(d, c));
          }
        } else {
          for (arma::uword d = 0; d < D; ++d) {
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(
                gamma_B(d, c).slice(s) * X_B(c).slice(i).col(0)
              );
            }
            B(d, c).each_slice() = Btmp.t();
            log_B(d, c) = arma::log(B(d, c));
          }
        }
      }
    }
  }
  void update_omega(const arma::uword i) {
    if (icpt_only_omega) {
      omega = softmax(gamma_omega.col(0));
    } else {
      omega = softmax(gamma_omega * X_omega.col(i));  
    }
    log_omega = arma::log(omega);
  }
  void update_log_py(const arma::uword i) {
    log_py.zeros();
    for (arma::uword d = 0; d < D; ++d) {
      for (arma::uword t = 0; t < Ti(i); ++t) {
        for (arma::uword c = 0; c < C; ++c) {
          log_py.slice(d).col(t) += log_B(d, c).slice(t).col(obs(c, t, i));
        }
      }
    }
  }
  void estep_pi(const arma::uword i, const arma::uword d, 
                const arma::vec& log_alpha, 
                const arma::vec& log_beta, const double ll) {
    E_pi(d).col(i) = arma::exp(log_alpha + log_beta - ll);
    // set minuscule values to zero in order to avoid numerical issues
    E_pi(d).col(i).clean(minval);
  }
  void estep_A(const arma::uword i, const arma::uword d, 
               const arma::mat& log_alpha, const arma::mat& log_beta, 
               const double ll) {
    for (arma::uword k = 0; k < S; ++k) { // from
      for (arma::uword j = 0; j < S; ++j) { // to
        for (arma::uword t = 0; t < (Ti(i) - 1); ++t) { // time
          E_A(d, k)(j, i, t + 1) = exp(log_alpha(k, t) + log_A(d)(k, j, t + 1) + 
            log_beta(j, t + 1) + log_py(j, t + 1, d) - ll);
        }
      }
      // set minuscule values to zero in order to avoid numerical issues
      E_A(d, k).col(i).clean(minval);
    }
  }
  void estep_B(const arma::uword i, const arma::uword d, 
               const arma::mat& log_alpha, const arma::mat& log_beta, 
               const double ll) {
    for (arma::uword k = 0; k < S; ++k) { // state
      for (arma::uword t = 0; t < Ti(i); ++t) { // time
        double pp = exp(log_alpha(k, t) + log_beta(k, t) - ll);
        for (arma::uword c = 0; c < C; ++c) { // channel
          if (obs(c, t, i) < M(c) && pp > minval) {
            E_B(d, c)(t, i, k) = pp;
          } else {
            E_B(d, c)(t, i, k) = 0.0;
          }
        }
      }
    }
  }
  void estep_omega(const arma::uword i, const arma::vec ll_i, 
                   const double ll) {
    E_omega.col(i) = arma::exp(ll_i - ll);
    // set minuscule values to zero in order to avoid numerical issues
    E_omega.col(i).clean(minval);
  }
  void mstep_pi(const arma::uword print_level);
  void mstep_A(const arma::uword print_level);
  void mstep_B(const arma::uword print_level);
  void mstep_omega(const arma::uword print_level);
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
  double objective_B(const arma::vec& x, arma::vec& grad);
  double objective_omega(const arma::vec& x, arma::vec& grad);
  static double objective_pi_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<mnhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_pi(x_vec, grad_vec);
  }
  static double objective_A_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<mnhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_A(x_vec, grad_vec);
  }
  static double objective_B_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<mnhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_B(x_vec, grad_vec);
  }
  static double objective_omega_wrapper(unsigned n, const double* x, double* grad, void* data) {
    auto* self = static_cast<mnhmm*>(data);
    arma::vec x_vec(const_cast<double*>(x), n, false, true);
    arma::vec grad_vec(grad, n, false, true);
    return self->objective_omega(x_vec, grad_vec);
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
        Rcpp::Named("eta_omega") = Rcpp::wrap(eta_omega),
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
    logp.fill(-arma::datum::inf);
    double logp_d;
    arma::uvec q_d(q.n_rows);
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_omega || i == 0) {
        update_omega(i);
      }
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
      for (arma::uword d = 0; d < D; ++d) {
        logp_d = univariate_viterbi(
          q_d,
          log_omega(d) + log_pi(d),
          log_A(d).slices(0, Ti(i) - 1),
          log_py.slice(d).cols(0, Ti(i) - 1)
        );
        if (logp_d > logp(i)) {
          logp(i) = logp_d;
          q.col(i) = q_d;
        }
      }
    }
  }
  void forward(arma::cube& log_alpha) {
    for (arma::uword i = 0; i < N; ++i) {
      if (!icpt_only_omega || i == 0) {
        update_omega(i);
      }
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
      for (arma::uword d = 0; d < D; ++d) {
        arma::subview<double> submat = 
          log_alpha.slice(i).rows(d * S, (d + 1) * S - 1);
        univariate_forward(
          submat,
          log_omega(d) + log_pi(d),
          log_A(d), 
          log_py.slice(d).cols(0, Ti(i) - 1)
        );
      }
    }
  }
  
  void backward(arma::cube& log_beta) {
    for (arma::uword i = 0; i < N; ++i) {
      if (iv_A || i == 0) {
        update_A(i);
      }
      if (arma::any(iv_B) || i == 0) {
        update_B(i);
      }
      update_log_py(i);
      for (arma::uword d = 0; d < D; ++d) {
        arma::subview<double> submat = 
          log_beta.slice(i).rows(d * S, (d + 1) * S - 1);
        univariate_backward(
          submat,
          log_A(d), 
          log_py.slice(d).cols(0, Ti(i) - 1)
        );
      }
    }
  }
};
#endif
