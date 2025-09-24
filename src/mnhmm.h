#ifndef MNHMM_H
#define MNHMM_H

#include "config.h"

class mnhmm { 
  
public:
  mnhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    arma::field<arma::mat>&& X_B,
    const arma::mat& X_omega,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool icpt_only_omega,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::field<arma::mat>& gamma_pi,
    const arma::field<arma::cube>& gamma_A,
    const arma::field<arma::cube>& gamma_B,
    const arma::mat& gamma_omega,
    double maxval_ = arma::datum::inf,
    double minval_ = -1.0);
  
  virtual ~mnhmm() = default;
  void update_pi(const arma::uword i);
  void update_A(const arma::uword i);
  virtual void update_B(const arma::uword i);
  void update_omega(const arma::uword i);
  void update_py(const arma::uword i);
  Rcpp::List viterbi();
  arma::vec loglik();
  arma::field<arma::mat> forward();
  arma::field<arma::mat> backward();
  Rcpp::List predict();
  Rcpp::List simulate();
  Rcpp::List log_objective(
      const arma::mat& Qs, 
      const arma::field<arma::mat>& Qm,
      const arma::mat& Qd
  );
  
  // data
  const arma::field<arma::umat>& obs;
  const arma::uvec& Ti;
  const arma::uvec& M;
  const arma::uword N;
  const arma::uword C;
  const arma::uword S;
  const arma::uword D;
  const arma::mat& X_pi;
  // field of length N containing K_A x T_i matrices
  const arma::field<arma::mat>& X_A;
  // N x C field containing K_B x T_i matrices
  arma::field<arma::mat> X_B;
  const arma::mat& X_omega;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const arma::uvec& icpt_only_B;
  const bool icpt_only_omega;
  const bool iv_A;
  const arma::uvec& iv_B;
  const bool tv_A; 
  const arma::uvec& tv_B;
  
  arma::field<arma::mat> gamma_pi;
  arma::field<arma::cube> gamma_A;
  arma::field<arma::cube> gamma_B;
  arma::mat gamma_omega;
  
  // pi, A, B, omega, and log_p(y) of _one_ id we are currently working with
  arma::cube py;
  arma::field<arma::vec> pi;
  arma::field<arma::cube> A;
  arma::field<arma::cube> B;
  arma::vec omega;
  arma::field<arma::vec> log_pi;
  arma::field<arma::cube> log_A;
  arma::field<arma::cube> log_B;
  arma::vec log_omega;
  
  const double maxval;
  const double minval;
  
private:
  void gradient_omega(
      arma::mat& grad, 
      const arma::vec& pcp, 
      const arma::uword i
  );
  void gradient_pi(
      arma::mat& grad, 
      arma::vec& tmpvec,
      const arma::vec& pcp,
      const arma::mat& beta, 
      const arma::uword i,
      const arma::uword d
  );
  void gradient_A(
      arma::mat& grad, 
      arma::vec& tmpvec1,
      arma::vec& tmpvec2,
      const arma::vec& pcp,
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword t, 
      const arma::uword s, 
      const arma::uword d
  );
  virtual void gradient_B_t1(
      arma::mat& grad, 
      arma::vec& tmpvec,
      const arma::vec& pcp,
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword s, 
      const arma::uword c, 
      const arma::uword d
  );
  void gradient_B(
      arma::mat& grad, 
      arma::vec& tmpvec, 
      const arma::vec& pcp,
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword s, 
      const arma::uword t, 
      const arma::uword c,
      const arma::uword d
  );
};
#endif
