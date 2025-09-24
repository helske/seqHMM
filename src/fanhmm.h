#ifndef FANHMM_H
#define FANHMM_H

#include "config.h"
#include "nhmm.h"

class fanhmm : public nhmm {
  
public:
  fanhmm(
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
    const arma::vec& prior_y,
    const Rcpp::List& W_X_B,
    double maxval_ = arma::datum::inf,
    double minval_ = -1.0);
  
  void update_B(const arma::uword i) override;
  Rcpp::List predict( 
      arma::field<arma::mat>&& W_A, 
      arma::field<arma::mat>&& W_B
  );
  Rcpp::List simulate(
      arma::field<arma::mat>&& W_A, 
      arma::field<arma::mat>&& W_B
  );
  const arma::vec& prior_y;
  const bool fixed_0;
  const arma::field<arma::vec> W_X_B;
  
private:
  void gradient_B_t1(
      arma::mat& grad,
      arma::vec& tmpvec,
      const arma::mat& beta,
      const arma::uword i,
      const arma::uword s,
      const arma::uword c
  ) override;
  
  arma::field<arma::cube> B1;
};

#endif
