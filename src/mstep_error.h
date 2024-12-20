#ifndef MSTEPERROR_H
#define MSTEPERROR_H

#include <RcppArmadillo.h>

template<typename Model>
Rcpp::List mstep_error_mnhmm(int return_code, const Model& model, 
                             int iter, double relative_change, 
                             double absolute_change, double absolute_x_change, 
                             double relative_x_change) {
  if (return_code != 0) {
    return Rcpp::List::create(
      Rcpp::Named("return_code") = return_code,
      Rcpp::Named("eta_omega") = Rcpp::wrap(model.eta_omega),
      Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
      Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
      Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
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

template<typename Model>
Rcpp::List mstep_error_nhmm(int return_code, const Model& model, 
                             int iter, double relative_change, 
                             double absolute_change, double absolute_x_change, 
                             double relative_x_change) {
  if (return_code != 0) {
    return Rcpp::List::create(
      Rcpp::Named("return_code") = return_code,
      Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
      Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
      Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
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

template<typename Model>
Rcpp::List mstep_error_fanhmm(int return_code, const Model& model, 
                            int iter, double relative_change, 
                            double absolute_change, double absolute_x_change, 
                            double relative_x_change) {
  if (return_code != 0) {
    return Rcpp::List::create(
      Rcpp::Named("return_code") = return_code,
      Rcpp::Named("eta_pi") = Rcpp::wrap(model.eta_pi),
      Rcpp::Named("eta_A") = Rcpp::wrap(model.eta_A),
      Rcpp::Named("eta_B") = Rcpp::wrap(model.eta_B),
      Rcpp::Named("rho_A") = Rcpp::wrap(model.rho_A),
      Rcpp::Named("rho_B") = Rcpp::wrap(model.rho_B),
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
#endif
