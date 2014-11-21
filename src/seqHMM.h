#ifndef SEQHMM_H
#define SEQHMM_H

#include <RcppArmadillo.h>

double logSumExp(const double& x,const double& y);

void internalForwardMC(const arma::mat& transition, const arma::cube& emission, 
                     const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha,arma::mat& scales);

void internalForward(const arma::mat& transition, const arma::mat& emission, 
                     const arma::vec& init, const arma::Mat<int>& obs, arma::cube& alpha,arma::mat& scales);

void internalBackwardMC(const arma::mat& transition, const arma::cube& emission, 
                      const arma::Cube<int>& obs, arma::cube& beta,const arma::mat& scales);
                       
void internalBackward(const arma::mat& transition, const arma::mat& emission, 
                      const arma::Mat<int>& obs, arma::cube& beta,const arma::mat& scales);
                      
                      
void internalForwardMCLog(const arma::mat& transition, const arma::cube& emission, 
                     const arma::vec& init, const arma::Cube<int>& obs, arma::cube& alpha);

void internalForwardLog(const arma::mat& transition, const arma::mat& emission, 
                     const arma::vec& init, const arma::Mat<int>& obs, arma::cube& alpha);

void internalBackwardMCLog(const arma::mat& transition, const arma::cube& emission, 
                      const arma::Cube<int>& obs, arma::cube& beta);
                       
void internalBackwardLog(const arma::mat& transition, const arma::mat& emission, 
                      const arma::Mat<int>& obs, arma::cube& beta);



#endif
