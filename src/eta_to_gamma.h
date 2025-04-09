#ifndef ETA_TO_GAMMA_H
#define ETA_TO_GAMMA_H

#include "config.h"

arma::mat eta_to_gamma(const arma::mat& eta, const arma::mat& Q);
arma::cube eta_to_gamma(const arma::cube& eta, const arma::mat& Q);
arma::field<arma::mat> eta_to_gamma(
    const arma::field<arma::mat>& eta, 
    const arma::mat& Q
);
arma::field<arma::cube> eta_to_gamma(
    const arma::field<arma::cube>& eta,
    const arma::mat& Q
);
arma::field<arma::cube> eta_to_gamma(
    const arma::field<arma::cube>& eta, 
    const arma::field<arma::mat>& Q, 
    const arma::uword D = 1
);

#endif


