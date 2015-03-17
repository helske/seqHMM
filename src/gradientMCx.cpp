//#include "seqHMM.h"
//using namespace Rcpp;
//
//// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::export]]
//
//NumericVector gradientMCx(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
//IntegerVector obsArray,NumericVector rowSumsA,NumericVector rowSumsB_,
//IntegerVector transNZ, IntegerVector emissNZ, NumericVector expPsi,NumericMatrix coef_, NumericMatrix X_) { 
//  
//
//  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
//  IntegerVector oDims = obsArray.attr("dim"); //k,n,r
//
//  arma::mat transition(transitionMatrix.begin(),eDims[0],eDims[0],false);
//  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1],eDims[2],false);
//  arma::icube obs(obsArray.begin(), oDims[0], oDims[1],oDims[2],false); 
//  arma::imat ANZ(transNZ.begin(),eDims[0],eDims[0],false);
//  arma::icube BNZ(emissNZ.begin(), eDims[0], eDims[1]-1,eDims[2],false);
//  arma::mat rowSumsB(rowSumsB_.begin(),eDims[0],eDims[2],false);
//  
//    int q = coef_.nrow();
//  arma::mat coef(coef_.begin(),q,eDims[0]); // qxm
//  arma::mat X(X_.begin(),oDims[0],q);       // kxq
//  arma::mat init = exp(X*coef).t(); // m x k
//  arma::rowvec sumInit = sum(init,0); // k
//  init.each_row() /= sumInit;
//  init = log(init); 
//  arma::cube alpha(eDims[0],oDims[1],oDims[0]); //m,n,k
//  arma::cube beta(eDims[0],oDims[1],oDims[0]); //m,n,k 
//  
//  arma::mat transitionLog = log(transition); 
//  arma::cube emissionLog = log(emission);
//  
//  internalForwardMCx(transitionLog, emissionLog, init, obs, alpha);
//  internalBackwardMC(transitionLog, emissionLog, obs, beta);     
//  
//  arma::vec ll(oDims[0]);
//  
//  double tmp=0.0;
//  double neginf = -arma::math::inf();
//  
//  for(int k=0;k<oDims[0];k++){
//    tmp =neginf;
//    for(int i = 0; i < eDims[0]; i++){
//      if(alpha(i,oDims[1]-1,k)>neginf){
//        tmp = logSumExp(alpha(i,oDims[1]-1,k),tmp); 
//      }
//    }
//    ll(k) = tmp;
//  }
//   arma::vec grad(expPsi.size()+q*(eDims[0]-1),arma::fill::zeros);
//  int countgrad=0;
//  // transitionMatrix
//  for(int i = 0; i < eDims[0]; i++){   
//    arma::uvec ind = arma::find(ANZ.row(i));
//    
//    if(ind.n_elem>1){ 
//      arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
//      for(unsigned int j = 0; j < ind.n_elem; j++){
//        for(int k = 0; k < oDims[0]; k++){ 
//          for(int t = 0; t < (oDims[1]-1); t++){  
//            tmp = 0.0;
//            for(int r=0; r < oDims[2]; r++){
//              tmp += emissionLog(ind(j),obs(k,t+1,r),r);
//            }
//            gradRow(j) += exp(alpha(i,t,k)+tmp+beta(ind(j),t+1,k)-ll(k));
//          }
//        }
//      }      
//      
//      for(unsigned int j = 0; j < ind.n_elem; j++){
//        if(ANZ(i,ind(j))!=2){
//          arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
//          dpsi(j) = 1.0;    
//          arma::uvec indi(1);
//          indi(0)=i;
//          dpsi = (dpsi-transition(indi,ind))*expPsi(countgrad)/rowSumsA(i);  
//          grad(countgrad) = arma::as_scalar(dpsi*gradRow);
//          countgrad ++;
//        }
//      }
//    }
//  }
//  
//  // emissionMatrix
//  for(int r=0; r < oDims[2]; r++){
//    for(int i = 0; i < eDims[0]; i++){   
//      arma::uvec ind = arma::find(BNZ.slice(r).row(i)); 
//      
//      if(ind.n_elem>1){ 
//        arma::vec gradRow(ind.n_elem,arma::fill::zeros);  
//        for(unsigned int j = 0; j < ind.n_elem; j++){
//          for(int k = 0; k < oDims[0]; k++){               
//            if(obs(k,0,r)==ind(j)){   
//               tmp = 0.0;
//              for(int r2=0; r2 < oDims[2]; r2++){
//                if(r2!=r){
//                  tmp += emissionLog(i,obs(k,0,r2),r2);
//                }
//              }
//              gradRow(j) += exp(init(i,k)+tmp+beta(i,0,k)-ll(k));
//            }       
//            for(int t = 0; t < (oDims[1]-1); t++){              
//              if(obs(k,t+1,r)==ind(j)){
//                tmp = 0.0;
//                for(int r2=0; r2 < oDims[2]; r2++){
//                  if(r2!=r){
//                    tmp += emissionLog(i,obs(k,t+1,r2),r2);
//                  }
//                }
//                gradRow(j) += arma::accu(exp(alpha.slice(k).col(t)+tmp+transitionLog.col(i)+beta(i,t+1,k)-ll(k)));
//              }          
//            }
//          }
//        }
//        for(unsigned int j = 0; j < ind.n_elem; j++){
//          if(BNZ(i,ind(j),r)!=2){
//            arma::rowvec dpsi(ind.n_elem,arma::fill::zeros);
//            dpsi(j) = 1.0;    
//            arma::uvec indi(1);
//            indi(0)=i;
//            dpsi = (dpsi-emission.slice(r).submat(indi,ind))*expPsi(countgrad)/rowSumsB(i,r);  
//            grad(countgrad) = arma::as_scalar(dpsi*gradRow);
//            countgrad ++;
//          }
//        }
//      }
//    }
//    
//  }
//   // InitProbs
//  for(int k = 0; k < oDims[0]; k++){ 
//    for(unsigned int j = 1; j < eDims[0]; j++){
//      double tmp = 0.0;
//      double tmp2 = 0.0;
//      for(int i = 0; i< eDims[0]; i++){
//        if(i!=j){
//          tmp += arma::as_scalar(exp(dot(coef.col(i),X.row(k))));
//          for(int r=0; r < oDims[2]; r++){
//          tmp2 += arma::as_scalar(exp(emissionLog(i,obs(k,0,r),r)+beta(i,0,k)-ll(k)));
//          }
//        }
//      }
//      for(int r=0; r < oDims[2]; r++){
//       tmp *= exp(emissionLog(j,obs(k,0,r),r)+beta(j,0,k)-ll(k));
//      }
//      grad.subvec(expPsi.size()+q*(j-1),expPsi.size()+q*j-1) += X.row(k).t()*
//      arma::as_scalar(exp(dot(coef.col(j),X.row(k)))/pow(sumInit(k),2)*
//      (tmp - tmp2));
//    }
//  }
//  
//  return wrap(grad);
//}