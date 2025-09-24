// EM algorithm for MNHMMs
#include "config.h"
#include "EM_mnhmm.h"
#include "create_Q.h"
#include "logsumexp.h"
#include "softmax.h"
#include "sum_to_zero.h"
#include "eta_to_gamma.h"
#include "forward.h"
#include "backward.h"
#include "mfanhmm.h"
#include <nloptrAPI.h>
#include <chrono>

EM_mnhmm::EM_mnhmm(
  mnhmm& model, 
  const arma::mat& Qs, 
  const arma::field<arma::mat>& Qm, 
  const arma::mat& Qd, 
  const double lambda,
  const arma::uword maxeval, 
  const double ftol_abs, 
  const double ftol_rel, 
  const double xtol_abs, 
  const double xtol_rel, 
  const arma::uword print_level,
  const arma::uword maxeval_m, 
  const double ftol_abs_m, 
  const double ftol_rel_m, 
  const double xtol_abs_m, 
  const double xtol_rel_m, 
  const arma::uword print_level_m,
  const double bound,
  const double tolg)
  : model(model), fan_model(dynamic_cast<const mfanhmm*>(&model)), 
    Qs(Qs), Qm(Qm), Qd(Qd), lambda(lambda),
    eta_pi(model.D), eta_A(model.D), eta_B(model.D, model.C), 
    eta_omega(Qd.t() * model.gamma_omega), E_pi(model.D),
    E_A(model.D, model.S), E_B(model.D, model.C), E_omega(model.D, model.N),
    opt_B(model.C, nullptr), 
    maxeval(maxeval), ftol_abs(ftol_abs), ftol_rel(ftol_rel), 
    xtol_abs(xtol_abs), xtol_rel(xtol_rel), print_level(print_level),
    maxeval_m(maxeval_m), ftol_abs_m(ftol_abs_m), ftol_rel_m(ftol_rel_m), 
    xtol_abs_m(xtol_abs_m), xtol_rel_m(xtol_rel_m), print_level_m(print_level_m),
    bound(bound), tolg(tolg)
{
  for (arma::uword d = 0; d < model.D; ++d) {
    eta_pi(d) = Qs.t() * model.gamma_pi(d);
    eta_A(d) = arma::cube(model.S - 1, model.X_A(0).n_rows, model.S);
    E_pi(d) = arma::mat(model.S, model.N);
    for (arma::uword s = 0; s < model.S; ++s) {
      eta_A(d).slice(s) = Qs.t() * model.gamma_A(d).slice(s);
      E_A(d, s) = arma::cube(model.S, model.N, model.Ti.max(), arma::fill::zeros);
    }
    for (arma::uword c = 0; c < model.C; ++c) {
      eta_B(d, c) = arma::cube(model.M(c) - 1, model.X_B(c, 0).n_rows, model.S);
      for (arma::uword s = 0; s < model.S; ++s) {
        eta_B(d, c).slice(s) = Qm(c).t() * model.gamma_B(d, c).slice(s);
      }
      E_B(d, c) = arma::cube(model.Ti.max(), model.N, model.S, arma::fill::zeros);
    }
  }
}
EM_mnhmm::~EM_mnhmm() {
  if (opt_pi) {
    nlopt_destroy(opt_pi);
  }
  if (opt_A) {
    nlopt_destroy(opt_A);
  }    
  for (auto& opt : opt_B) {
    if (opt) {
      nlopt_destroy(opt);
      opt = nullptr;
    }
  }
  if (opt_omega) {
    nlopt_destroy(opt_omega);
  }
}
void EM_mnhmm::update_gamma_pi() {
  model.gamma_pi = eta_to_gamma(eta_pi, Qs);
}
void EM_mnhmm::update_gamma_A() {
  model.gamma_A = eta_to_gamma(eta_A, Qs);
}
void EM_mnhmm::update_gamma_B() {
  model.gamma_B = eta_to_gamma(eta_B, Qm, model.D);
}
void EM_mnhmm::update_gamma_omega() {
  model.gamma_omega = eta_to_gamma(eta_omega, Qd);
}

Rcpp::List EM_mnhmm::mstep_error(int iter, 
                                 double relative_change, 
                                 double absolute_change, 
                                 double absolute_x_change, 
                                 double relative_x_change) {
  
  if (mstep_return_code != 0) {
    return Rcpp::List::create(
      Rcpp::Named("return_code") = mstep_return_code,
      Rcpp::Named("eta_pi") = Rcpp::wrap(eta_pi),
      Rcpp::Named("eta_A") = Rcpp::wrap(eta_A),
      Rcpp::Named("eta_B") = Rcpp::wrap(eta_B),
      Rcpp::Named("eta_omega") = Rcpp::wrap(eta_omega),
      Rcpp::Named("gamma_pi") = Rcpp::wrap(model.gamma_pi),
      Rcpp::Named("gamma_A") = Rcpp::wrap(model.gamma_A),
      Rcpp::Named("gamma_B") = Rcpp::wrap(model.gamma_B),
      Rcpp::Named("gamma_omega") = Rcpp::wrap(model.gamma_omega),
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

void EM_mnhmm::estep_pi(const arma::uword i, const arma::uword d, 
                        const arma::vec& alpha, const arma::vec& beta, 
                        const double pcp, const double scale1) {
  E_pi(d).col(i) = pcp * alpha % beta / scale1;
}

void EM_mnhmm::estep_A(const arma::uword i, const arma::uword d, 
                       const arma::mat& alpha, const arma::mat& beta, 
                       const double pcp) {
  
  for (arma::uword k = 0; k < model.S; ++k) { // from
    for (arma::uword j = 0; j < model.S; ++j) { // to
      for (arma::uword t = 0; t < (model.Ti(i) - 1); ++t) { // time
        E_A(d, k)(j, i, t + 1) = pcp * alpha(k, t) * model.A(d)(k, j, t + 1) * 
          beta(j, t + 1) * model.py(j, t + 1, d);
      }
    }
  }
}

void EM_mnhmm::estep_B(const arma::uword i, const arma::uword d, 
                       const arma::mat& alpha, const arma::mat& beta, 
                       const double pcp, const arma::vec& scales) {
  
  for (arma::uword k = 0; k < model.S; ++k) { // state
    for (arma::uword t = 0; t < model.Ti(i); ++t) { // time
      double pp = alpha(k, t) * beta(k, t) / scales(t);
      for (arma::uword c = 0; c < model.C; ++c) { // channel
        if (model.obs(i)(c, t) < model.M(c) && pp > model.minval) {
          E_B(d, c)(t, i, k) = pcp * pp;
        } else {
          E_B(d, c)(t, i, k) = 0.0;
        }
      }
    }
  }
}
void EM_mnhmm::estep_omega(const arma::uword i, const arma::vec& likelihood) {
  
  E_omega.col(i) = likelihood;
}

void EM_mnhmm::mstep_pi() {
  
  mstep_return_code = 0;
  // Use closed form solution
  if (model.icpt_only_pi && lambda < 1e-12) {
    arma::vec tmp(model.S);
    for (arma::uword d = 0; d < model.D; ++d) {
      tmp = arma::sum(E_pi(d), 1);
      tmp.clamp(std::numeric_limits<double>::min(), arma::datum::inf);
      eta_pi(d) = Qs.t() * log(tmp);
      if (!eta_pi(d).is_finite()) {
        mstep_return_code = -100;
        return;
      }
    }
    return;
  }
  if (!opt_pi) {
    Rcpp::stop("Optimizer opt_pi not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_pi, objective_pi_wrapper, this);
  arma::vec x(eta_pi(0).memptr(), eta_pi.n_elem, false, true);
  double minf;
  int return_code;
  arma::vec grad(eta_pi(0).n_elem);
  for (arma::uword d = 0; d < model.D; ++d) {
    current_d = d;
    arma::vec x(eta_pi(d).memptr(), eta_pi(d).n_elem, false, true);
    double ll = objective_pi(x, grad);
    last_val = std::numeric_limits<double>::infinity();
    abs_change = 0;
    rel_change = 0;
    mstep_iter = 0;
    if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
      return_code = 1; // already converged (L-BFGS gradient tolerance)
    } else {
      return_code = nlopt_optimize(opt_pi, x.memptr(), &minf);
      // nlopt_optimize can return generic failure code, check that we still have progressed
      if (return_code == -1 && (abs_change < ftol_abs || rel_change < ftol_rel)) {
        return_code = 7;
      }
    }
    if (print_level_m > 0) {
      Rcpp::Rcout<<"M-step of initial probabilities in cluster "<<d + 1
                 <<" ended with return code "<<return_code
                 <<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
      if (print_level_m > 1) {
        Rcpp::Rcout<<"Relative change "<<rel_change<<", absolute change "<<abs_change<<std::endl;
      }
    }
    if (return_code < 0) {
      mstep_return_code = return_code - 110;
      return;
    }
  }
}
void EM_mnhmm::mstep_A() {
  
  mstep_return_code = 0;
  // Use closed form solution
  if (model.icpt_only_A && lambda < 1e-12) {
    arma::vec tmp(model.S);
    for (arma::uword d = 0; d < model.D; ++d) {
      for (arma::uword s = 0; s < model.S; ++s) { // from
        for (arma::uword k = 0; k < model.S; ++k) { //to
          tmp(k) = arma::accu(E_A(d, s).row(k));
        }
        tmp.clamp(std::numeric_limits<double>::min(), arma::datum::inf);
        eta_A(d).slice(s).col(0) = Qs.t() * log(tmp);
        if (!eta_A(d).slice(s).col(0).is_finite()) {
          mstep_return_code = -200;
          return;
        }
      }
    }
    return;
  }
  if (!opt_A) {
    Rcpp::stop("Optimizer opt_A not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_A, objective_A_wrapper, this);
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_A(0).slice(0).n_elem);
  for (arma::uword d = 0; d < model.D; ++d) {
    current_d = d;
    for (arma::uword s = 0; s < model.S; ++s) {
      current_s = s;
      arma::vec x(eta_A(d).slice(s).memptr(), eta_A(d).slice(s).n_elem, false, true);
      ll = objective_A(x, grad);
      last_val = std::numeric_limits<double>::infinity();
      abs_change = 0;
      rel_change = 0;
      mstep_iter = 0;
      if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
        return_code = 1; // already converged (L-BFGS gradient tolerance)
      } else {
        return_code = nlopt_optimize(opt_A, x.memptr(), &minf);
        // nlopt_optimize can return generic failure code, check that we still have progressed
        if (return_code == -1 && (abs_change < ftol_abs || rel_change < ftol_rel)) {
          return_code = 7;
        }
      }
      if (print_level_m > 0) {
        Rcpp::Rcout<<"M-step of transition probabilities of state "<<s + 1
                   <<" in cluster "<<d + 1
                   <<" ended with return code "<< return_code
                   <<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
        if (print_level_m > 1) {
          Rcpp::Rcout<<"Relative change "<<rel_change<<", absolute change "<<abs_change<<std::endl;
        }
      }
      if (return_code < 0) {
        mstep_return_code = return_code - 210;
        return;
      }
    }
  }
}
void EM_mnhmm::mstep_B() {
  mstep_return_code = 0;
  // use closed form solution
  if (arma::all(model.icpt_only_B) && lambda < 1e-12) {
    for (arma::uword c = 0; c < model.C; ++c) {
      arma::vec tmp(model.M(c));
      for (arma::uword d = 0; d < model.D; ++d) {
        for (arma::uword s = 0; s < model.S; ++s) {
          tmp.zeros();
          for (arma::uword i = 0; i < model.N; ++i) {
            for (arma::uword t = 0; t < model.Ti(i); ++t) {
              if (model.obs(i)(c, t) < model.M(c)) {
                tmp(model.obs(i)(c, t)) += E_B(d, c)(t, i, s);
              }
            }
          }
          tmp.clamp(std::numeric_limits<double>::min(), arma::datum::inf);
          eta_B(d, c).slice(s).col(0) = Qm(c).t() * log(tmp);
          if (!eta_B(d, c).slice(s).col(0).is_finite()) {
            mstep_return_code = -300;
            return;
          }
        }
      }
    }
    return;
  }
  double minf;
  int return_code;
  double ll;
  for (arma::uword c = 0; c < model.C; ++c) {
    current_c = c;
    arma::vec grad(eta_B(0, c).slice(0).n_elem);
    nlopt_set_min_objective(opt_B[c], objective_B_wrapper, this);
    for (arma::uword d = 0; d < model.D; ++d) {
      current_d = d;
      for (arma::uword s = 0; s < model.S; ++s) {
        current_s = s;
        arma::vec x(eta_B(d, c).slice(s).memptr(), eta_B(d, c).slice(s).n_elem, false, true);
        ll = objective_B(x, grad);
        last_val = std::numeric_limits<double>::infinity();
        abs_change = 0;
        rel_change = 0;
        mstep_iter = 0;
        if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
          return_code = 1; // already converged (L-BFGS gradient tolerance)
        } else {
          return_code = nlopt_optimize(opt_B[c], x.memptr(), &minf);
          // nlopt_optimize can return generic failure code, check that we still have progressed
          if (return_code == -1 && (abs_change < ftol_abs || rel_change < ftol_rel)) {
            return_code = 7;
          }
        }
        if (print_level_m > 0) {
          Rcpp::Rcout<<"M-step of emission probabilities of state "<<s + 1
                     <<" in cluster "<<d + 1
                     <<" of response "<<c + 1
                     <<" ended with return code "<< return_code
                     <<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
          if (print_level_m > 1) {
            Rcpp::Rcout<<"Relative change "<<rel_change<<", absolute change "<<abs_change<<std::endl;
          }
        }
        if (return_code < 0) {
          mstep_return_code = return_code - 310;
          return;
        }
      }
    }
  }
}
void EM_mnhmm::mstep_omega() {
  mstep_return_code = 0;
  // Use closed form solution
  if (model.icpt_only_omega && lambda < 1e-12) {
    arma::vec tmp = arma::sum(E_omega, 1);
    tmp.clamp(std::numeric_limits<double>::min(), arma::datum::inf);
    eta_omega = Qd.t() * log(tmp);
    if (!eta_omega.is_finite()) {
      mstep_return_code = -400;
      return;
    }
    return;
  }
  if (!opt_omega) {
    Rcpp::stop("Optimizer opt_omega not initialized! Shouldn't be possible, file an issue.");
  }
  nlopt_set_min_objective(opt_omega, objective_omega_wrapper, this);
  arma::vec x(eta_omega.memptr(), eta_omega.n_elem, false, true);
  
  double minf;
  int return_code;
  double ll;
  arma::vec grad(eta_omega.n_elem);
  ll = objective_omega(x, grad);
  last_val = std::numeric_limits<double>::infinity();
  abs_change = 0;
  rel_change = 0;
  mstep_iter = 0;
  if (arma::norm(grad, "inf") < 1e-8 && std::isfinite(ll)) {
    return_code = 1; // already converged (L-BFGS gradient tolerance)
  } else {
    return_code = nlopt_optimize(opt_omega, x.memptr(), &minf);
    // nlopt_optimize can return generic failure code, check that we still have progressed
    if (return_code == -1 && (abs_change < ftol_abs || rel_change < ftol_rel)) {
      return_code = 7;
    }
  }
  if (print_level_m > 0) {
    Rcpp::Rcout<<"M-step of cluster probabilities ended with return code "<<
      return_code<<" after "<<mstep_iter + 1<<" iterations."<<std::endl;
    if (print_level_m > 1) {
      Rcpp::Rcout<<"Relative change "<<rel_change<<", absolute change "<<abs_change<<std::endl;
    }
  }
  if (return_code < 0) {
    mstep_return_code = return_code - 410;
  }
}

double EM_mnhmm::objective_pi_wrapper(unsigned n, const double* x, double* grad, void* data) {
  auto* self = static_cast<EM_mnhmm*>(data);
  arma::vec x_vec(const_cast<double*>(x), n, false, true);
  arma::vec grad_vec(grad, n, false, true);
  return self->objective_pi(x_vec, grad_vec);
}

double EM_mnhmm::objective_A_wrapper(unsigned n, const double* x, double* grad, void* data) {
  auto* self = static_cast<EM_mnhmm*>(data);
  arma::vec x_vec(const_cast<double*>(x), n, false, true);
  arma::vec grad_vec(grad, n, false, true);
  return self->objective_A(x_vec, grad_vec);
}

double EM_mnhmm::objective_B_wrapper(unsigned n, const double* x, double* grad, void* data) {
  auto* self = static_cast<EM_mnhmm*>(data);
  arma::vec x_vec(const_cast<double*>(x), n, false, true);
  arma::vec grad_vec(grad, n, false, true);
  return self->objective_B(x_vec, grad_vec);
}

double EM_mnhmm::objective_omega_wrapper(unsigned n, const double* x, double* grad, void* data) {
  auto* self = static_cast<EM_mnhmm*>(data);
  arma::vec x_vec(const_cast<double*>(x), n, false, true);
  arma::vec grad_vec(grad, n, false, true);
  return self->objective_omega(x_vec, grad_vec);
}

double EM_mnhmm::objective_pi(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  eta_pi(current_d) = arma::mat(x.memptr(), model.S - 1, model.X_pi.n_rows);
  model.gamma_pi(current_d) = sum_to_zero(eta_pi(current_d), Qs);
  grad.zeros();
  arma::mat tQs = Qs.t();
  for (arma::uword i = 0; i < model.N; ++i) {
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    const arma::vec& counts = E_pi(current_d).col(i);
    double sum_epi = arma::accu(counts);
    double val = arma::dot(counts, model.log_pi(current_d));
    value -= val;
    grad -= arma::vectorise(tQs * (counts - sum_epi * model.pi(current_d)) * model.X_pi.col(i).t());
  }
  value += 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
  value /= model.N;
  grad += lambda * x;
  grad /= model.N;
  abs_change = value - last_val;
  rel_change = std::abs(abs_change) / (std::abs(last_val) + 1e-12);
  last_val = value;
  return value;
}

double EM_mnhmm::objective_A(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::mat eta_Arow = arma::mat(x.memptr(), model.S - 1, model.X_A(0).n_rows);
  arma::mat gamma_Arow = sum_to_zero(eta_Arow, Qs);
  arma::vec A1(model.S);
  arma::vec log_A1(model.S);
  grad.zeros();
  if (!model.iv_A && !model.tv_A) {
    A1 = softmax(gamma_Arow * model.X_A(0).col(0));
    log_A1 = log(A1);
  }
  arma::mat tQs = Qs.t();
  for (arma::uword i = 0; i < model.N; ++i) {
    if (model.iv_A && !model.tv_A) {
      A1 = softmax(gamma_Arow * model.X_A(i).col(0));
      log_A1 = log(A1);
    }
    for (arma::uword t = 1; t < model.Ti(i); ++t) {
      const arma::vec& counts = E_A(current_d, current_s).slice(t).col(i);
      if (model.tv_A) {
        A1 = softmax(gamma_Arow * model.X_A(i).col(t));
        log_A1 = log(A1);
      }
      double sum_ea = arma::accu(counts);
      double val = arma::dot(counts, log_A1);
      value -= val;
      grad -= arma::vectorise(tQs * (counts - sum_ea * A1) * model.X_A(i).col(t).t());
    }
  }
  value += 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
  value /= model.N;
  grad += lambda * x;
  grad /= model.N;
  abs_change = value - last_val;
  rel_change = std::abs(abs_change) / (std::abs(last_val) + 1e-12);
  last_val = value;
  return value;
}

double EM_mnhmm::objective_B(const arma::vec& x, arma::vec& grad) {
  mstep_iter++;
  double value = 0;
  arma::uword Mc = model.M(current_c);
  arma::mat eta_Brow = arma::mat(x.memptr(), Mc - 1, model.X_B(current_c, 0).n_rows);
  arma::mat gamma_Brow = sum_to_zero(eta_Brow, Qm(current_c));
  arma::vec B1(Mc);
  arma::vec log_B1(Mc);
  grad.zeros();
  if (!model.iv_B(current_c) && !model.tv_B(current_c)) {
    B1 = softmax(gamma_Brow * model.X_B(current_c, 0).col(0));
    log_B1 = log(B1);
  }
  arma::mat I(Mc, Mc, arma::fill::eye);
  arma::mat tQm = Qm(current_c).t();
  
  // Check if model is a fanhmm
  arma::uword J = 0;
  bool marginalize_y0 = false;
  if (fan_model) {
    J = fan_model->prior_y.n_elem;
    if (J > 1) {
      marginalize_y0 = true;
    }
  }
  arma::mat B1j(Mc, J);
  arma::uword y;
  for (arma::uword i = 0; i < model.N; ++i) {
    if (!marginalize_y0 && model.iv_B(current_c) && !model.tv_B(current_c)) {
      B1 = softmax(gamma_Brow * model.X_B(current_c, i).col(0));
      log_B1 = log(B1);
    }
    if (marginalize_y0) {
      y = model.obs(i)(current_c, 0);
      if (y < Mc) {
        double e_b = E_B(current_d, current_c)(0, i, current_s);
        if (e_b > 0) {
          B1.zeros();
          for (arma::uword j = 0; j < J; ++j) {
            B1j.col(j) = softmax(gamma_Brow * fan_model->W_X_B(j, current_c, i));
            B1 += B1j.col(j) * fan_model->prior_y(j);
          }
          log_B1 = log(B1);
          double val = e_b * log_B1(y);
          value -= val;
          for (arma::uword j = 0; j < J; ++j) {
            grad -= arma::vectorise(tQm * e_b  * B1j(y, j) / B1(y) * (I.col(y) - B1j.col(j)) *
              fan_model->W_X_B(j, current_c, i).t() * fan_model->prior_y(j));
          }
        }
      }
    }
    for (arma::uword t = 0 + marginalize_y0; t < model.Ti(i); ++t) {
      arma::uword y = model.obs(i)(current_c, t);
      if (y < Mc) {
        double e_b = E_B(current_d, current_c)(t, i, current_s);
        if (e_b > 0) {
          if (model.tv_B(current_c)) {
            B1 = softmax(gamma_Brow * model.X_B(current_c, i).col(t));
            log_B1 = log(B1);
          }
          double val = e_b * log_B1(y);
          value -= val;
          grad -= arma::vectorise(tQm * e_b  * (I.col(y) - B1) * 
            model.X_B(current_c, i).col(t).t());
        }
      }
    }
  }
  value += 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
  value /= model.N;
  grad += lambda * x;
  grad /= model.N;
  abs_change = value - last_val;
  rel_change = std::abs(abs_change) / (std::abs(last_val) + 1e-12);
  last_val = value;
  return value;
}

double EM_mnhmm::objective_omega(const arma::vec& x, arma::vec& grad) {
  
  mstep_iter++;
  double value = 0;
  eta_omega = arma::mat(x.memptr(), model.D - 1, model.X_omega.n_rows);
  model.gamma_omega = sum_to_zero(eta_omega, Qd);
  grad.zeros();
  arma::mat tQd = Qd.t();
  for (arma::uword i = 0; i < model.N; ++i) {
    if (!model.icpt_only_omega || i == 0) {
      model.update_omega(i);
    }
    const arma::vec& counts = E_omega.col(i);
    double val = arma::dot(counts, model.log_omega);
    value -= val;
    grad -= arma::vectorise(tQd * (counts - model.omega) * model.X_omega.col(i).t());
  }
  value += 0.5 * lambda * std::pow(arma::norm(x, 2), 2);
  value /= model.N;
  grad += lambda * x;
  grad /= model.N;
  abs_change = value - last_val;
  rel_change = std::abs(abs_change) / (std::abs(last_val) + 1e-12);
  last_val = value;
  return value;
}

Rcpp::List EM_mnhmm::run() {
  
  const arma::uword n_obs = arma::accu(model.Ti);
  arma::uword n_omega = eta_omega.n_elem;
  arma::uword n_pi = eta_pi(0).n_elem;
  arma::uword n_A = eta_A(0).n_elem;
  arma::uvec n_Bc(model.C);
  for (arma::uword c = 0; c < model.C; ++c) {
    n_Bc(c) = eta_B(0, c).n_elem;
  }
  arma::uword n_B = arma::accu(n_Bc);
  opt_pi = nlopt_create(NLOPT_LD_LBFGS, eta_pi(0).n_elem);
  nlopt_set_xtol_abs1(opt_pi, xtol_abs_m);
  nlopt_set_ftol_abs(opt_pi, ftol_abs_m);
  nlopt_set_xtol_rel(opt_pi, xtol_rel_m);
  nlopt_set_ftol_rel(opt_pi, ftol_rel_m);
  nlopt_set_maxeval(opt_pi, maxeval_m);
  nlopt_set_lower_bounds1(opt_pi, -bound);
  nlopt_set_upper_bounds1(opt_pi, bound);
  nlopt_set_vector_storage(opt_pi, 10);
  // nlopt_set_param(opt_pi, "tolg", tolg);
  
  opt_A = nlopt_create(NLOPT_LD_LBFGS, eta_A(0).slice(0).n_elem);
  nlopt_set_xtol_abs1(opt_A, xtol_abs_m);
  nlopt_set_ftol_abs(opt_A, ftol_abs_m);
  nlopt_set_xtol_rel(opt_A, xtol_rel_m);
  nlopt_set_ftol_rel(opt_A, ftol_rel_m);
  nlopt_set_maxeval(opt_A, maxeval_m);
  nlopt_set_lower_bounds1(opt_A, -bound);
  nlopt_set_upper_bounds1(opt_A, bound);
  nlopt_set_vector_storage(opt_A, 10);
  // nlopt_set_param(opt_A, "tolg", tolg);
  
  for (arma::uword c = 0; c < model.C; ++c) {
    opt_B[c] = nlopt_create(NLOPT_LD_LBFGS, eta_B(0, c).slice(0).n_elem);
    nlopt_set_xtol_abs1(opt_B[c], xtol_abs_m);
    nlopt_set_ftol_abs(opt_B[c], ftol_abs_m);
    nlopt_set_xtol_rel(opt_B[c], xtol_rel_m);
    nlopt_set_ftol_rel(opt_B[c], ftol_rel_m);
    nlopt_set_maxeval(opt_B[c], maxeval_m);
    nlopt_set_lower_bounds1(opt_B[c], -bound);
    nlopt_set_upper_bounds1(opt_B[c], bound);
    nlopt_set_vector_storage(opt_B[c], 10);
    // nlopt_set_param(opt_B[c], "tolg", tolg);
  }
  opt_omega = nlopt_create(NLOPT_LD_LBFGS, eta_omega.n_elem);
  nlopt_set_xtol_abs1(opt_omega, xtol_abs_m);
  nlopt_set_ftol_abs(opt_omega, ftol_abs_m);
  nlopt_set_xtol_rel(opt_omega, xtol_rel_m);
  nlopt_set_ftol_rel(opt_omega, ftol_rel_m);
  nlopt_set_maxeval(opt_omega, maxeval_m);
  nlopt_set_lower_bounds1(opt_omega, -bound);
  nlopt_set_upper_bounds1(opt_omega, bound);
  nlopt_set_vector_storage(opt_omega, 10);
  // nlopt_set_param(opt_omega, "tolg", tolg);
  
  arma::rowvec new_pars(n_omega + model.D * (n_pi + n_A + n_B));
  arma::rowvec pars(n_omega + model.D * (n_pi + n_A + n_B));
  
  // EM algorithm begins
  pars.cols(0, n_omega - 1) = arma::vectorise(eta_omega).t();
  arma::uword ii = n_omega;
  for (arma::uword d = 0; d < model.D; ++d) {
    pars.cols(ii, ii + n_pi - 1) = arma::vectorise(eta_pi(d)).t();
    ii += n_pi;
  }
  for (arma::uword d = 0; d < model.D; ++d) {
    pars.cols(ii, ii + n_A - 1) = arma::vectorise(eta_A(d)).t();
    ii += n_A;
  }
  for (arma::uword c = 0; c < model.C; ++c) {
    for (arma::uword d = 0; d < model.D; ++d) {
      pars.cols(ii, ii + n_Bc(c) - 1) = arma::vectorise(eta_B(d, c)).t();
      ii += n_Bc(c);
    }
  }
  double relative_change = ftol_rel + 1.0;
  double absolute_change = ftol_abs + 1.0;
  double relative_x_change = xtol_rel + 1.0;
  double absolute_x_change = xtol_abs+ 1.0;
  arma::uword iter = 0;
  double ll_new;
  double ll;
  arma::cube alpha(model.S, model.Ti.max(), model.D);
  arma::cube beta(model.S, model.Ti.max(), model.D);
  arma::mat scales(model.Ti.max(), model.D);
  arma::vec loglik(model.N);
  arma::vec loglik_i(model.D);
  arma::vec lls(model.D);
  // Initial log-likelihood
  for (arma::uword i = 0; i < model.N; ++i) {
    if (!model.icpt_only_pi || i == 0) {
      model.update_pi(i);
    }
    if (model.iv_A || i == 0) {
      model.update_A(i);
    }
    if (arma::any(model.iv_B) || i == 0) {
      model.update_B(i);
    }
    if (!model.icpt_only_omega || i == 0) {
      model.update_omega(i);
    }
    model.update_py(i);
    for (arma::uword d = 0; d < model.D; ++d) {
      arma::subview_col<double> scales_col = scales.col(d);
      loglik_i(d) = univariate_forward(
        alpha.slice(d), scales_col, model.pi(d), model.A(d), 
        model.py.slice(d), model.Ti(i)
      );
      univariate_backward(
        beta.slice(d), scales.col(d), model.A(d), 
        model.py.slice(d), model.Ti(i)
      );
    }
    lls = model.log_omega + loglik_i;
    loglik(i) = logSumExp(lls);
    lls = arma::exp(lls - loglik(i));
    for (arma::uword d = 0; d < model.D; ++d) {
      estep_pi(i, d, alpha.slice(d).col(0), beta.slice(d).col(0), lls(d), scales(0, d));
      estep_A(i, d, alpha.slice(d), beta.slice(d), lls(d));
      estep_B(i, d, alpha.slice(d), beta.slice(d), lls(d), scales.col(d));
    }
    estep_omega(i, lls);
  }
  ll = arma::accu(loglik);
  double penalty_term = 0.5 * lambda * std::pow(arma::norm(pars, 2), 2);
  ll -= penalty_term;
  ll /= n_obs;
  if (print_level > 0) {
    Rcpp::Rcout<<"Initial value of the log-likelihood: "<<ll<<std::endl;
    if (print_level > 1) {
      Rcpp::Rcout<<"Initial parameter values"<<std::endl;
      Rcpp::Rcout<<pars<<std::endl;
    }
  }
  // check for user interrupt every two seconds
  auto start_time = std::chrono::steady_clock::now();
  const std::chrono::seconds check_interval(2);
  
  while (relative_change > ftol_rel && absolute_change > ftol_abs &&
         absolute_x_change > xtol_abs && 
         relative_x_change > xtol_rel && iter < maxeval) {
    
    auto current_time = std::chrono::steady_clock::now();
    if (current_time - start_time >= check_interval) {
      Rcpp::checkUserInterrupt();
      start_time = current_time; // Reset the timer
    }
    
    iter++;
    ll_new = 0;
    mstep_pi();
    if (mstep_return_code != 0) {
      return mstep_error(
        iter, relative_change, absolute_change, absolute_x_change, 
        relative_x_change
      );
    }  
    mstep_A();
    if (mstep_return_code != 0) {
      return mstep_error(
        iter, relative_change, absolute_change, absolute_x_change, 
        relative_x_change
      );
    }
    mstep_B();
    if (mstep_return_code != 0) {
      return mstep_error(
        iter, relative_change, absolute_change, absolute_x_change, 
        relative_x_change
      );
    }
    mstep_omega();
    if (mstep_return_code != 0) {
      return mstep_error(
        iter, relative_change, absolute_change, absolute_x_change, 
        relative_x_change
      );
    }
    // Update model
    update_gamma_pi();
    update_gamma_A();
    update_gamma_B();
    update_gamma_omega();
    for (arma::uword i = 0; i < model.N; ++i) {
      if (!model.icpt_only_pi || i == 0) {
        model.update_pi(i);
      }
      if (model.iv_A || i == 0) {
        model.update_A(i);
      }
      if (arma::any(model.iv_B) || i == 0) {
        model.update_B(i);
      }
      if (!model.icpt_only_omega || i == 0) {
        model.update_omega(i);
      }
      model.update_py(i);
      for (arma::uword d = 0; d < model.D; ++d) {
        arma::subview_col<double> scales_col = scales.col(d);
        loglik_i(d) = univariate_forward(
          alpha.slice(d), scales_col, model.pi(d), model.A(d), model.py.slice(d), 
          model.Ti(i)
        );
        univariate_backward(
          beta.slice(d), scales.col(d), model.A(d), 
          model.py.slice(d), model.Ti(i)
        );
      }
      lls = model.log_omega + loglik_i;
      loglik(i) = logSumExp(lls);
      lls = arma::exp(lls - loglik(i));
      for (arma::uword d = 0; d < model.D; ++d) {
        estep_pi(i, d, alpha.slice(d).col(0), beta.slice(d).col(0), lls(d), scales(0, d));
        estep_A(i, d, alpha.slice(d), beta.slice(d), lls(d));
        estep_B(i, d, alpha.slice(d), beta.slice(d), lls(d), scales.col(d));
      }
      estep_omega(i, lls);
    }
    ll_new = arma::accu(loglik);
    new_pars.cols(0, n_omega - 1) = arma::vectorise(eta_omega).t();
    ii = n_omega;
    for (arma::uword d = 0; d < model.D; ++d) {
      new_pars.cols(ii, ii + n_pi - 1) = arma::vectorise(eta_pi(d)).t();
      ii += n_pi;
    }
    for (arma::uword d = 0; d < model.D; ++d) {
      new_pars.cols(ii, ii + n_A - 1) = arma::vectorise(eta_A(d)).t();
      ii += n_A;
    }
    for (arma::uword c = 0; c < model.C; ++c) {
      for (arma::uword d = 0; d < model.D; ++d) {
        new_pars.cols(ii, ii + n_Bc(c) - 1) = arma::vectorise(eta_B(d, c)).t();
        ii += n_Bc(c);
      }
    }
    penalty_term = 0.5 * lambda * std::pow(arma::norm(new_pars, 2), 2);
    ll_new -= penalty_term;
    ll_new /= n_obs;
    relative_change = (ll_new - ll) / (std::abs(ll) + 1e-12);
    absolute_change = ll_new - ll;
    absolute_x_change = arma::max(arma::abs(new_pars - pars));
    relative_x_change = arma::norm(new_pars - pars, 1) / arma::norm(pars, 1);
    if (print_level > 0) {
      Rcpp::Rcout<<"Iteration: "<<iter<<std::endl;
      Rcpp::Rcout<<"           "<<"log-likelihood: "<<ll_new<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of log-likelihood: "<<relative_change<<std::endl;
      Rcpp::Rcout<<"           "<<"absolute change of scaled log-likelihood: "<<absolute_change<<std::endl;
      Rcpp::Rcout<<"           "<<"relative change of parameters: "<<relative_x_change<<std::endl;
      Rcpp::Rcout<<"           "<<"maximum absolute change of parameters: "<<absolute_x_change<<std::endl;
      if (print_level > 1) {
        Rcpp::Rcout << "current parameter values"<< std::endl;
        Rcpp::Rcout<<new_pars<<std::endl;
      }
    }
    ll = ll_new;
    pars = new_pars;
    if (absolute_change < -1e-8) {
      Rcpp::warning("EM algorithm encountered decreasing log-likelihood.");
    }
  }
  int return_code = 0;
  if (iter >= maxeval) {
    return_code = 5;
  } else if (relative_change < ftol_rel || absolute_change < ftol_abs) {
    return_code = 3;
  } else if (relative_x_change < xtol_rel || absolute_x_change < xtol_abs) {
    return_code = 4;
  }
  update_gamma_pi();
  update_gamma_A();
  update_gamma_B();
  update_gamma_omega();
  return Rcpp::List::create(
    Rcpp::Named("return_code") = return_code,
    Rcpp::Named("eta_pi") = Rcpp::wrap(eta_pi),
    Rcpp::Named("eta_A") = Rcpp::wrap(eta_A),
    Rcpp::Named("eta_B") = Rcpp::wrap(eta_B),
    Rcpp::Named("eta_omega") = Rcpp::wrap(eta_omega),
    Rcpp::Named("gamma_pi") = Rcpp::wrap(model.gamma_pi),
    Rcpp::Named("gamma_A") = Rcpp::wrap(model.gamma_A),
    Rcpp::Named("gamma_B") = Rcpp::wrap(model.gamma_B),
    Rcpp::Named("gamma_omega") = Rcpp::wrap(model.gamma_omega),
    Rcpp::Named("logLik") = ll * n_obs,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("relative_f_change") = relative_change,
    Rcpp::Named("absolute_f_change") = absolute_change,
    Rcpp::Named("absolute_x_change") = absolute_x_change,
    Rcpp::Named("relative_x_change") = relative_x_change
  );
}
