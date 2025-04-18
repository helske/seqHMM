#ifndef FANHMM_H
#define FANHMM_H

#include "config.h"
#include "nhmm.h"

class fanhmm : public nhmm {
  
public:
  fanhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
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
      const arma::field<arma::cube>& W_A, 
      const arma::field<arma::cube>& W_B
  );
  Rcpp::List simulate(
      const arma::field<arma::cube>& W_A, 
      const arma::field<arma::cube>& W_B
  );
  const arma::vec& prior_y;
  const bool fixed_0;
  const arma::field<arma::cube> W_X_B;
};

#endif
