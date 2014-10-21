
#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logLikMCHMM(NumericVector transitionMatrix, NumericVector emissionArray, 
NumericVector initialProbs, IntegerVector obsArray) {  
  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  
  arma::colvec init(initialProbs.begin(),eDims[0]);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0]);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2]);
  arma::Cube<int> obs(obsArray.begin(), oDims[0], oDims[1], oDims[2]);
  
  
  transition = log(transition); 
  emission = log(emission); 
  init = log(init); 
  
  arma::vec alpha(eDims[0]); //m,n,k
  arma::vec alphatmp(eDims[0]); //m,n,k
  //arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  
  
  double ll=0.0;
  
  double sumtmp;
  double tmp;
  double neginf = -arma::math::inf();
  
//  for(int k = 0; k < oDims[0]; k++){
//    for(int i=0; i < eDims[0]; i++){      
//      alpha.at(i,0,k) = init.at(i);
//      for(int r = 0; r < oDims[2]; r++){
//        alpha.at(i,0,k) += emission.at(i,obs.at(k,0,r),r);
//      }
//    }
//  }
//  
//  
//  for(int k = 0; k < oDims[0]; k++){
//    for(int t = 1; t < oDims[1]; t++){  
//      for(int i = 0; i < eDims[0]; i++){
//        sumtmp = neginf;
//        for(int j = 0; j < eDims[0]; j++){
//          tmp = alpha.at(j,t-1,k) + transition.at(j,i);
//          if(tmp > neginf){
//            sumtmp = logSumExp(sumtmp,tmp);
//          }
//        }
//        
//        for(int r = 0; r < oDims[2]; r++){
//          sumtmp += emission.at(i,obs.at(k,t,r),r);
//        }
//        alpha.at(i,t,k) = sumtmp;
//      }
//      
//    }
//  }
//  
//
//  
//  arma::vec llk(oDims[0]);
//  for(int k=0;k<oDims[0];k++){
//    tmp =neginf;
//    for(int i = 0; i < eDims[0]; i++){
//      if(alpha(i,oDims[1]-1,k)>neginf){
//        tmp = logSumExp(alpha(i,oDims[1]-1,k),tmp); 
//      }
//    }
//    llk(k) = tmp;// + log(sum(exp(alpha.slice(k).col(oDims[1]-1)-tmp)));
//  }
//  
  
  //  for(int k=0;k<oDims[0];k++){   
  //    tmp = max(alpha.slice(k).col(oDims[1]-1)); 
  //    ll(k) = tmp + log(sum(exp(alpha.slice(k).col(oDims[1]-1)-tmp)));
  //  }
  
  //ll = accu(llk);
  
    
    for(int k = 0; k < oDims[0]; k++){    
      
      
      for(int i=0; i < eDims[0]; i++){      
        alpha(i) = init(i);
        for(int r = 0; r < oDims[2]; r++){
          alpha(i) += emission(i,obs(k,0,r),r);
        }
      }    
     
      for(int t = 1; t < oDims[1]; t++){  
        for(int i = 0; i < eDims[0]; i++){
          sumtmp = neginf;
          for(int j = 0; j < eDims[0]; j++){
            tmp = alpha(j) + transition(j,i);
            if(tmp > neginf){
              sumtmp = logSumExp(sumtmp,tmp);
            }
          }
          
          for(int r = 0; r < oDims[2]; r++){
            sumtmp += emission(i,obs(k,t,r),r);
          }
          alphatmp(i) = sumtmp;
        }
        alpha = alphatmp;
       
      }
      
      
      tmp = neginf;
      for(int i = 0; i < eDims[0]; i++){
        if(alpha(i)>neginf){
          tmp = logSumExp(alpha(i),tmp); 
        }
      }
      ll += tmp;    
    
    }
    
  return ll;
}

