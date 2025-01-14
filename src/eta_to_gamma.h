#ifndef ETATOGAMMA_H
#define ETATOGAMMA_H

#include <RcppArmadillo.h>
arma::mat eta_to_gamma(const arma::mat& eta);
arma::cube eta_to_gamma(const arma::cube& eta);
arma::field<arma::mat> eta_to_gamma(const arma::field<arma::mat>& eta);
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta);
arma::mat eta_to_gamma(const arma::mat& eta, const arma::mat& Q);
arma::cube eta_to_gamma(const arma::cube& eta, const arma::mat& Q);
arma::field<arma::mat> eta_to_gamma(const arma::field<arma::mat>& eta, const arma::mat& Q);
arma::field<arma::cube> eta_to_gamma(const arma::field<arma::cube>& eta, const arma::mat& Q);
#endif
