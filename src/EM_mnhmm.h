#ifndef EM_MNHMM_H
#define EM_MNHMM_H

#include "mnhmm.h"
#include "mfanhmm.h"
#include "create_Q.h"
#include <nloptrAPI.h>

class EM_mnhmm {
  
public:
  EM_mnhmm(mnhmm& model, 
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
           const double tolg
  );
  ~EM_mnhmm();
  Rcpp::List run();
  
private:
  // functions
  
  void update_gamma_pi();
  void update_gamma_A();
  void update_gamma_B();
  void update_gamma_omega();
  
  Rcpp::List mstep_error(
      int iter, 
      double relative_change, 
      double absolute_change, 
      double absolute_x_change, 
      double relative_x_change
  );
  
  void estep_pi(
      const arma::uword i, 
      const arma::uword d, 
      const arma::vec& alpha, 
      const arma::vec& beta, 
      const double pcp,
      const double scale
  );
  void estep_A(
      const arma::uword i, 
      const arma::uword d, 
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const double pcp
  );
  void estep_B(
      const arma::uword i, 
      const arma::uword d, 
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const double pcp,
      const arma::vec& scales
  );
  
  void estep_omega(const arma::uword i, const arma::vec& likelihood);
  void mstep_pi();
  void mstep_A();
  void mstep_B();
  void mstep_omega();
  
  double objective_pi(const arma::vec& x, arma::vec& grad);
  double objective_A(const arma::vec& x, arma::vec& grad);
  double objective_B(const arma::vec& x, arma::vec& grad);
  double objective_omega(const arma::vec& x, arma::vec& grad);
  
  static double objective_pi_wrapper(unsigned n, const double* x, double* grad, void* data);
  static double objective_A_wrapper(unsigned n, const double* x, double* grad, void* data);
  static double objective_B_wrapper(unsigned n, const double* x, double* grad, void* data);
  static double objective_omega_wrapper(unsigned n, const double* x, double* grad, void* data);
  
  // data
  mnhmm& model;
  const mfanhmm* fan_model;
  const arma::mat Qs;
  const arma::field<arma::mat> Qm;
  const arma::mat Qd;
  const double lambda;
  
  // coefficients //
  arma::field<arma::mat> eta_pi;
  arma::field<arma::cube> eta_A;
  arma::field<arma::cube> eta_B;  
  arma::mat eta_omega;
  
  // excepted counts
  arma::field<arma::mat> E_pi;
  arma::field<arma::cube> E_A;
  arma::field<arma::cube> E_B;
  arma::mat E_omega;
  
  arma::uword current_s = 0;
  arma::uword current_c = 0; 
  arma::uword current_d = 0;
  unsigned int mstep_iter = 0;
  int mstep_return_code = 0;  
  nlopt_opt opt_pi = nullptr;
  nlopt_opt opt_A = nullptr;
  std::vector<nlopt_opt> opt_B;
  nlopt_opt opt_omega = nullptr;
  
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
