#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector gradientMCx(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB_, NumericVector sumInit,
IntegerVector transNZ, IntegerVector emissNZ, IntegerVector initNZ, NumericVector expPsi,
NumericMatrix coefs, NumericMatrix X_, IntegerVector numberOfStates) { 

  
  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
  
  arma::vec init(initialProbs.begin(),eDims[0],false);
  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1],eDims[2],false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1],oDims[2],false); 
  arma::imat ANZ(transNZ.begin(),eDims[0],eDims[0],false);
  arma::icube BNZ(emissNZ.begin(), eDims[0], eDims[1]-1,eDims[2],false);
  arma::ivec INZ(initNZ.begin(), eDims[0],false);
  arma::mat rowSumsB(rowSumsB_.begin(),eDims[0],eDims[2],false);
  
  int q = coefs.nrow();
  arma::mat coef(coefs.begin(),q,numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(),oDims[0],q);
  
  arma::mat lweights = exp(X*coef).t();
  if(!lweights.is_finite()){
    return -std::numeric_limits<double>::max();
  }
  arma::rowvec sumweights = sum(lweights,0);
  
  lweights.each_row() /= sumweights;
  lweights = log(lweights); 
  
  arma::vec initLog = log(init); 
  
  
  arma::mat initk(eDims[0],oDims[0]);
  
  for(int k = 0; k < oDims[0]; k++){    
    initk.col(k) = initLog + reparma(lweights.col(k),numberOfStates);
  }
  
  
  arma::mat transitionLog = log(transition); 
  arma::cube emissionLog = log(emission);
  
  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k 
  
  internalForwardMCx(transitionLog, emissionLog, initk, obs, alpha);
  internalBackwardMC(transitionLog, emissionLog, obs, beta);     
  
  arma::vec ll(oDims[0]);
  
  double tmp=0.0;
  double neginf = -arma::math::inf();
  
  for(int k=0;k<oDims[0];k++){
    tmp =neginf;
    for(int i = 0; i < eDims[0]; i++){
      if(alpha(i,oDims[1]-1,k)>neginf){
        tmp = logSumExp(alpha(i,oDims[1]-1,k),tmp); 
      }
    }
    ll(k) = tmp;
  }
  arma::vec grad(expPsi.size()+q*(numberOfStates.size()-1),arma::fill::zeros);
  int countgrad=0;
  // transitionMatrix
  for(int i = 0; i < eDims[0]; i++){   
    arma::uvec ind = arma::find(ANZ.row(i));
    
    if(ind.n_elem>1){ 
      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
      for(unsigned int j = 0; j < ind.n_elem; j++){
        for(int k = 0; k < oDims[0]; k++){ 
          for(int t = 0; t < (oDims[1]-1); t++){  
            tmp = 0.0;
            for(int r=0; r < oDims[2]; r++){
              tmp += emissionLog(ind(j),obs(k,t+1,r),r);
            }
            gradRow(j) += exp(alpha(i,t,k)+tmp+beta(ind(j),t+1,k)-ll(k));
          }
        }
      }      
      
      for(unsigned int j = 0; j < ind.n_elem; j++){
        if(ANZ(i,ind(j))!=2){
          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
          dpsi(j) = 1.0;    
          arma::uvec indi(1);
          indi(0)=i;
          dpsi = (dpsi-transition(indi,ind))*expPsi(countgrad)/rowSumsA(i);  
          grad(countgrad) = arma::as_scalar(dpsi*gradRow);
          countgrad ++;
        }
      }
    }
  }
  
  // emissionMatrix
  for(int r=0; r < oDims[2]; r++){
    for(int i = 0; i < eDims[0]; i++){   
      arma::uvec ind = arma::find(BNZ.slice(r).row(i)); 
      
      if(ind.n_elem>1){ 
        arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
        for(unsigned int j = 0; j < ind.n_elem; j++){
          for(int k = 0; k < oDims[0]; k++){               
            if(obs(k,0,r)==ind(j)){   
              tmp = 0.0;
              for(int r2=0; r2 < oDims[2]; r2++){
                if(r2!=r){
                  tmp += emissionLog(i,obs(k,0,r2),r2);
                }
              }
              gradRow(j) += exp(initk(i,k)+tmp+beta(i,0,k)-ll(k));
            }       
            for(int t = 0; t < (oDims[1]-1); t++){              
              if(obs(k,t+1,r)==ind(j)){
                tmp = 0.0;
                for(int r2=0; r2 < oDims[2]; r2++){
                  if(r2!=r){
                    tmp += emissionLog(i,obs(k,t+1,r2),r2);
                  }
                }
                gradRow(j) += arma::accu(exp(alpha.slice(k).col(t)+tmp+transitionLog.col(i)+beta(i,t+1,k)-ll(k)));
              }          
            }
          }
        }
        for(unsigned int j = 0; j < ind.n_elem; j++){
          if(BNZ(i,ind(j),r)!=2){
            arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
            dpsi(j) = 1.0;    
            arma::uvec indi(1);
            indi(0)=i;
            dpsi = (dpsi-emission.slice(r).submat(indi,ind))*expPsi(countgrad)/rowSumsB(i,r);  
            grad(countgrad) = arma::as_scalar(dpsi*gradRow);
            countgrad ++;
          }
        }
      }
    }
    
  }
  // InitProbs  
 
        IntegerVector cumsumstate = cumsum(numberOfStates);
      for(unsigned int jj = 0; jj < numberOfStates.size(); jj++){
        arma::uvec ind = arma::find(INZ.subvec(cumsumstate(jj)-numberOfStates(jj),
        cumsumstate(jj)-1));  
        if(ind.n_elem>1){     
          arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
          for(unsigned int j = 0; j < ind.n_elem; j++){        
            for(int k = 0; k < oDims[0]; k++){ 
              tmp = 0.0;
              for(int r=0; r < oDims[2]; r++){
                tmp += emissionLog(cumsumstate(jj)-numberOfStates(jj)+j,obs(k,0,r),r);
              }
              gradRow(j) += exp(tmp+beta(cumsumstate(jj)-numberOfStates(jj)+j,0,k)-ll(k)+lweights(jj,k));
            }
          }
          
          for(unsigned int j = 0; j < ind.n_elem; j++){
            if(INZ(cumsumstate(jj)-numberOfStates(jj)+ind(j))!=2){
              arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
              dpsi(j) = 1.0;    
              dpsi = (dpsi-init(cumsumstate(jj)-numberOfStates(jj)+ind).t())*expPsi(countgrad)/sumInit(jj);
              grad(countgrad) = arma::as_scalar(dpsi*gradRow);          
              countgrad ++;
              
            }
          }
        }
      }
  //    arma::uvec ind = arma::find(INZ);    
  //  if(ind.n_elem>1){ 
  //    arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
  //    for(unsigned int j = 0; j < ind.n_elem; j++){
  //      for(int k = 0; k < oDims[0]; k++){ 
  //        tmp = 0.0;
  //        for(int r=0; r < oDims[2]; r++){
  //          tmp += emissionLog(j,obs(k,0,r),r);
  //        }
  //        gradRow(j) += exp(tmp+beta(j,0,k)-ll(k));
  //      }
  //    }
  //    for(unsigned int j = 0; j < ind.n_elem; j++){
  //      if(INZ(ind(j))!=2){
  //        arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
  //        dpsi(j) = 1.0;    
  //        dpsi = (dpsi-init(ind).t())*expPsi(countgrad)/sum(sumInit);  
  //        grad(countgrad) = arma::as_scalar(dpsi*gradRow);
  //        countgrad ++;
  //      }
  //    }
  //  }
    //  IntegerVector cumsumstate = cumsum(numberOfStates);
    //  for(unsigned int jj = 0; jj < numberOfStates.size(); jj++){
    //    arma::uvec ind = arma::find(INZ.subvec(cumsumstate(jj)-numberOfStates(0),
    //    cumsumstate(jj)-numberOfStates(0)+numberOfStates(jj)-1));  
    //    if(ind.n_elem>1){     
    //      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
    //      for(unsigned int j = 0; j < ind.n_elem; j++){        
    //        for(int k = 0; k < oDims[0]; k++){ 
    //          tmp = 0.0;
    //          for(int r=0; r < oDims[2]; r++){
    //            tmp += emissionLog(cumsumstate(jj)-numberOfStates(0)+j,obs(k,0,r),r);
    //          }
    //          gradRow(j) += exp(tmp+beta(cumsumstate(jj)-numberOfStates(0)+j,0,k)-ll(k));
    //        }
    //      }
    //      
    //      for(unsigned int j = 0; j < ind.n_elem; j++){
    //        if(INZ(cumsumstate(jj)-numberOfStates(0)+ind(j))!=2){
    //          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
    //          dpsi(j) = 1.0;    
    //          dpsi = (dpsi-init(cumsumstate(jj)-numberOfStates(0)+ind).t())*expPsi(countgrad)/sum(sumInit); //0 is wrong  
    //          grad(countgrad) = arma::as_scalar(dpsi*gradRow);          
    //          countgrad ++;
    //          
    //        }
    //      }
    //    }
    //  }
    // beta
    
    double tmp2;
    for(unsigned int jj = 1; jj < numberOfStates.size(); jj++){
      for(int k = 0; k < oDims[0]; k++){
        tmp2 = 0.0;
        for(unsigned int j = 0; j< numberOfStates.size(); j++){
          if(j!=jj){
            tmp2 += exp(dot(coef.col(j),X.row(k)));
          }
        }
        for(unsigned int j = 0; j < eDims[0]; j++){                
          tmp = 0.0;
          for(int r=0; r < oDims[2]; r++){
            tmp += emissionLog(j,obs(k,0,r),r);
          }        
          if(j>=(cumsumstate(jj)-numberOfStates(jj)) & j<cumsumstate(jj)){
            grad.subvec(expPsi.size()+q*(jj-1),expPsi.size()+q*jj-1) += 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*tmp2/sumweights(k);
          } else {
            grad.subvec(expPsi.size()+q*(jj-1),expPsi.size()+q*jj-1) -= 
            exp(tmp+beta(j,0,k)-ll(k)+initk(j,k))*X.row(k).t()*
            exp(dot(coef.col(jj),X.row(k)))/sumweights(k);
          }
        }
      }
    }
    //    arma::uvec ind = arma::find(INZ.subvec(cumsumstate(jj)-numberOfStates(0),
    //    cumsumstate(jj)-numberOfStates(0)+numberOfStates(jj)-1));  
    //    if(ind.n_elem>1){     
    //      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
    //      for(unsigned int j = 0; j < ind.n_elem; j++){        
    //        for(int k = 0; k < oDims[0]; k++){ 
    //          tmp = 0.0;
    //          for(int r=0; r < oDims[2]; r++){
    //            tmp += emissionLog(cumsumstate(jj)-numberOfStates(0)+j,obs(k,0,r),r);
    //          }
    //          gradRow(j) += exp(tmp+beta(cumsumstate(jj)-numberOfStates(0)+j,0,k)-ll(k));
    //        }
    //      }
    //      
    //      for(unsigned int j = 0; j < ind.n_elem; j++){
    //        if(INZ(cumsumstate(jj)-numberOfStates(0)+ind(j))!=2){
    //          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
    //          dpsi(j) = 1.0;    
    //          dpsi = (dpsi-init(cumsumstate(jj)-numberOfStates(0)+ind).t())*expPsi(countgrad)/sum(sumInit); //0 is wrong  
    //          grad(countgrad) += arma::as_scalar(dpsi*gradRow)*;                   
    //          
    //        }
    //      }
    //    }
    //     countgrad ++;
    //  }
    //arma::mat llc(oDims[0],numberOfStates.size());
    
    //  for(int k=0;k<oDims[0];k++){
    //    // all groups
    //    for(unsigned int jj = 0; jj < numberOfStates.size(); jj++){
    //      // log(P(seq_i|model_jj))
    //      tmp =neginf;
    //      for(int i = cumsumstate(jj)-numberOfStates(0); i < cumsumstate(jj)-numberOfStates(0)+numberOfStates(jj); i++){
    //        if(alpha(i,oDims[1]-1,k)>neginf){
    //          tmp = logSumExp(alpha(i,oDims[1]-1,k),tmp); 
    //        }
    //      }
    //      llc(k,jj) = tmp;
    //    }
    // all sequences
    //  for(int k = 0; k < oDims[0]; k++){ 
    //    // all groups
    //    for(unsigned int jj = 0; jj < numberOfStates.size(); jj++){ 
    //      
    //      double tmp = 0.0;
    //      double tmp2 = 0.0;
    //      for(int i = 0; i< numberOfStates.size(); i++){
    //        if(i!=j){
    //          tmp += arma::as_scalar(exp(dot(coef.col(i),X.row(k))));
    //          for(int r=0; r < oDims[2]; r++){
    //            tmp2 += arma::as_scalar(exp(emissionLog(i,obs(k,0,r),r)+beta(i,0,k)-llc(k,jj)));
    //          }
    //        }
    //      }
    //      for(int r=0; r < oDims[2]; r++){
    //        tmp *= exp(emissionLog(j,obs(k,0,r),r)+beta(j,0,k)-ll(k));
    //      }
    //      grad.subvec(expPsi.size()+q*(j-1),expPsi.size()+q*j-1) += X.row(k).t()*
    //      arma::as_scalar(exp(dot(coef.col(j),X.row(k)))/pow(sumInit,2)*(tmp - tmp2));
    //    }
    //  }
    
    return wrap(grad);
  }