#include "seqHMM.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// install_github( "Rcpp11/attributes" ) ; require('attributes') 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List hardEM(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
  IntegerVector obsArray, IntegerVector nSymbols, int itermax=100, double tol=1e-8, int trace=0) {  
  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(), eDims[0], true);
  arma::mat transition(transitionMatrix.begin(), eDims[0], eDims[0], true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  
  transition = log(transition); 
  emission = log(emission);
  init = log(init); 
  
  arma::umat q(oDims[0], oDims[1]);
  arma::vec logp(oDims[0]);
  
  viterbiForEM(transition, emission, init, obs, logp, q);
  
  double sumlogp = sum(logp);
  
  if(trace>0){
    Rcout<<"Log-likelihood of initial model: "<< sumlogp<<std::endl;
  }
  //  
  //  //EM-algorithm begins
  //  
  double change = tol+1.0;
  int iter = 0;
  double tmp;
  arma::mat ksii(eDims[0],eDims[0]);
  arma::cube gamma(eDims[0],eDims[1],eDims[2]);
  arma::vec delta(eDims[0]);
  
  while((change>tol) & (iter<itermax)){   
    iter++;
    gamma.zeros();
    ksii.zeros();
    delta.zeros();
    
    for(int k = 0; k < oDims[0]; k++){
      
      delta(q(k,0))++;     
      
      for(int t=0; t < (oDims[1]-1); t++){
        ksii(q(k,t),q(k,t+1))++;
        for(int r = 0; r < eDims[2]; r++){
          gamma(q(k,t),obs(k,t,r),r)++;  
        }
        
      }
      for(int r = 0; r < eDims[2]; r++){
        gamma(q(k,oDims[1]-1),obs(k,oDims[1]-1,r),r)++;  
      }
    }
    delta /= arma::as_scalar(arma::accu(delta));
    ksii.each_col() /= sum(ksii,1);  
    
    init = log(delta);
    transition = log(ksii);
    for(int r=0; r < eDims[2]; r++){
      gamma.slice(r).cols(0,nSymbols(r)-1).each_col() /= sum(gamma.slice(r).cols(0,nSymbols(r)-1),1);
      emission.slice(r).cols(0,nSymbols(r)-1) = log(gamma.slice(r).cols(0,nSymbols(r)-1));
    }
    
    
    viterbiForEM(transition, emission, init, obs, logp, q);
    tmp = sum(logp);
    change = (tmp - sumlogp)/(abs(sumlogp)+0.1);
    sumlogp = tmp;
    if(trace>1){
      Rcout<<"iter: "<< iter;
      Rcout<<" logLik: "<< sumlogp;
      Rcout<<" relative change: "<<change<<std::endl;
    }
    
  }
  if(trace>0){
    if(iter==itermax){
      Rcpp::Rcout<<"EM algorithm stopped after reaching the maximum number of "<<iter<<" iterations."<<std::endl;     
    } else{
      Rcpp::Rcout<<"EM algorithm stopped after reaching the relative change of "<<change;
      Rcpp::Rcout<<" after "<<iter<<" iterations."<<std::endl;
    }
    Rcpp::Rcout<<"Final log-likelihood: "<< sumlogp<<std::endl;
  }
  return List::create(Named("initialProbs") = wrap(exp(init)), Named("transitionMatrix") = wrap(exp(transition)),
    Named("emissionArray") = wrap(exp(emission)),Named("logLik") = sumlogp,Named("iterations")=iter,Named("change")=change);
}
