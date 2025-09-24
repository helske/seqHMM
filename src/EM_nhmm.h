#ifndef EM_NHMM_H
#define EM_NHMM_H

#include "nhmm.h"
#include "fanhmm.h"
#include <nloptrAPI.h>

class EM_nhmm {
public:
  EM_nhmm(
    nhmm& model, 
    const arma::mat& Qs, 
    const arma::field<arma::mat>& Qm,
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
    const double tolg);
  ~EM_nhmm();
  Rcpp::List run();
  
private:
  // functions
  void update_gamma_pi();
  void update_gamma_A();
  void update_gamma_B();
  
  Rcpp::List mstep_error(
      int iter, 
      double relative_change, 
      double absolute_change, 
      double absolute_x_change, 
      double relative_x_change
  );
  
  void estep_pi(
      const arma::uword i, 
      const arma::vec& alpha, 
      const arma::vec& beta,
      const double scale
  );
  void estep_A(
      const arma::uword i, 
      const arma::mat& alpha, 
      const arma::mat& beta
  );
  void estep_B(
      const arma::uword i, 
      const arma::mat& alpha, 
      const arma::mat& beta,
      const arma::vec& scales
  );
  
  void mstep_pi();
  void mstep_A();
  void mstep_B();
  
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
  double objective_B(const arma::vec& x, arma::vec& grad);
  
  static double objective_pi_wrapper(
      unsigned n, const double* x, double* grad, void* data)
    ;
  static double objective_A_wrapper(
      unsigned n, const double* x, double* grad, void* data
  );
  static double objective_B_wrapper(
      unsigned n, const double* x, double* grad, void* data
  );
  
  // data
  nhmm& model;
  const fanhmm* fan_model;
  const arma::mat Qs;
  const arma::field<arma::mat> Qm;
  const double lambda;
  
  //coefficients //
  arma::mat eta_pi;
  arma::cube eta_A;
  arma::field<arma::cube> eta_B;
  
  // excepted counts
  arma::mat E_pi;
  arma::field<arma::cube> E_A;
  arma::field<arma::cube> E_B;
  arma::uword current_s = 0;
  arma::uword current_c = 0; 
  unsigned int mstep_iter = 0;
  int mstep_return_code = 0;
  nlopt_opt opt_pi = nullptr;
  nlopt_opt opt_A = nullptr;
  std::vector<nlopt_opt> opt_B;
  
  const arma::uword maxeval; 
  const double ftol_abs; 
  const double ftol_rel; 
  const double xtol_abs; 
  const double xtol_rel; 
  const arma::uword print_level;
  const arma::uword maxeval_m; 
  const double ftol_abs_m; 
  const double ftol_rel_m; 
  const double xtol_abs_m; 
  const double xtol_rel_m; 
  const arma::uword print_level_m;
  const double bound;
  const double tolg;
  
  double last_val = std::numeric_limits<double>::infinity();
  double abs_change = 0;
  double rel_change = 0;
};

#endif
