#ifndef NHMM_H
#define NHMM_H

#include "config.h"

class nhmm {
  
public:
  
  nhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    arma::field<arma::mat>&& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& gamma_pi,
    const arma::cube& gamma_A,
    const arma::field<arma::cube>& gamma_B,
    double maxval_ = arma::datum::inf,
    double minval_ = -1.0);
  
  virtual ~nhmm() = default;
  
  void update_pi(const arma::uword i);
  void update_A(const arma::uword i);
  virtual void update_B(const arma::uword i);
  void update_py(const arma::uword i);
  Rcpp::List viterbi();
  arma::vec loglik();
  arma::field<arma::mat> forward();
  arma::field<arma::mat> backward();
  Rcpp::List predict();
  Rcpp::List simulate();
  Rcpp::List log_objective(
      const arma::mat& Qs, 
      const arma::field<arma::mat>& Qm
  );
  
  // data
  // field of length N containing C x T_i matrices
  const arma::field<arma::umat>& obs;
  const arma::uvec& Ti;
  const arma::uvec& M;
  const arma::uword N;
  const arma::uword C;
  const arma::uword S;
  const arma::mat& X_pi;
  // field of length N containing K_A x T_i matrices
  const arma::field<arma::mat>& X_A;
  // N x C field containing K_B x T_i matrices
  arma::field<arma::mat> X_B;
  const bool icpt_only_pi;
  const bool icpt_only_A;
  const arma::uvec& icpt_only_B;
  const bool iv_A;
  const arma::uvec& iv_B;
  const bool tv_A; 
  const arma::uvec& tv_B;
  
  arma::mat gamma_pi;
  arma::cube gamma_A;
  arma::field<arma::cube> gamma_B;
  
  // pi, A, B, and p(y) of _one_ id we are currently working with
  arma::mat py;
  arma::vec pi;
  arma::cube A;
  arma::field<arma::cube> B;
  arma::vec log_pi;
  arma::cube log_A;
  arma::field<arma::cube> log_B;
  const double maxval;
  const double minval;
  
private:
  
  void gradient_pi(
      arma::mat& grad, 
      arma::vec& tmpvec, 
      const arma::mat& beta, 
      const arma::uword i
  );
  void gradient_A(
      arma::mat& grad,
      arma::vec& tmpvec1,
      arma::vec& tmpvec2,
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword t, 
      const arma::uword s
  );
  virtual void gradient_B_t1(
      arma::mat& grad, 
      arma::vec& tmpvec, 
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword s, 
      const arma::uword c
  );
  void gradient_B(
      arma::mat& grad, 
      arma::vec& tmpvec, 
      const arma::mat& alpha, 
      const arma::mat& beta, 
      const arma::uword i, 
      const arma::uword s, 
      const arma::uword t, 
      const arma::uword c
  );
};

#endif
