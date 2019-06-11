// log-likelihood and gradients of HMM using log-space

#include "log_forward_backward.h"
#include "logsumexp.h"

// [[Rcpp::export]]
Rcpp::List log_objective(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, const arma::umat& ANZ,
  const arma::ucube& BNZ, const arma::uvec& INZ, arma::uvec& nSymbols, unsigned int threads) {
  
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), arma::fill::zeros);
  
  arma::vec initLog = log(init);
  arma::mat transitionLog = log(transition);
  arma::cube emissionLog = log(emission);
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  
  log_internalForward(transitionLog, emissionLog, initLog, obs, alpha, threads);
  log_internalBackward(transitionLog, emissionLog, obs, beta, threads);
  
  arma::vec ll(obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
  }
  
  arma::mat gradmat(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), obs.n_slices,
    arma::fill::zeros);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) default(shared)
    for (unsigned int k = 0; k < obs.n_slices; k++) {
      int countgrad = 0;
      
      // transitionMatrix
      arma::vec gradArow(emission.n_rows);
      arma::mat gradA(emission.n_rows, emission.n_rows);
      
      for (unsigned int i = 0; i < emission.n_rows; i++) {
        arma::uvec ind = arma::find(ANZ.row(i));
        
        if (ind.n_elem > 0) {
          gradArow.zeros();
          gradA.eye();
          gradA.each_row() -= transition.row(i);
          gradA.each_col() %= transition.row(i).t();
          
          for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
            for (unsigned int j = 0; j < emission.n_rows; j++) {
              double tmp = 0.0;
              for (unsigned int r = 0; r < obs.n_rows; r++) {
                tmp += emissionLog(j, obs(r, t + 1, k), r);
              }
              gradArow(j) += exp(alpha(i, t, k) + tmp + beta(j, t + 1, k) - ll(k));
            }
            
          }
          
          gradArow = gradA * gradArow;
          gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradArow.rows(ind);
          countgrad += ind.n_elem;
        }
      }
      // emissionMatrix
      for (unsigned int r = 0; r < obs.n_rows; r++) {
        arma::vec gradBrow(nSymbols[r]);
        arma::mat gradB(nSymbols[r], nSymbols[r]);
        for (unsigned int i = 0; i < emission.n_rows; i++) {
          arma::uvec ind = arma::find(BNZ.slice(r).row(i));
          if (ind.n_elem > 0) {
            gradBrow.zeros();
            gradB.eye();
            gradB.each_row() -= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1);
            gradB.each_col() %= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1).t();
            for (unsigned int j = 0; j < nSymbols[r]; j++) {
              if (obs(r, 0, k) == j) {
                double tmp = 0.0;
                for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                  if (r2 != r) {
                    tmp += emissionLog(i, obs(r2, 0, k), r2);
                  }
                }
                gradBrow(j) += exp(initLog(i) + tmp + beta(i, 0, k) - ll(k));
              }
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                if (obs(r, t + 1, k) == j) {
                  double tmp = 0.0;
                  for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                    if (r2 != r) {
                      tmp += emissionLog(i, obs(r2, t + 1, k), r2);
                    }
                  }
                  gradBrow(j) += arma::accu(
                    exp(alpha.slice(k).col(t) + tmp + transitionLog.col(i) 
                    + beta(i, t + 1, k) - ll(k)));
                }
              }
              
            }
            gradBrow = gradB * gradBrow;
            gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradBrow.rows(ind);
            countgrad += ind.n_elem;
            
          }
        }
      }
      // InitProbs
      arma::uvec ind = arma::find(INZ);
      if (ind.n_elem > 0) {
        arma::vec gradIrow(emission.n_rows);
        arma::mat gradI(emission.n_rows, emission.n_rows);
        
        gradIrow.zeros();
        gradI.zeros();
        gradI.eye();
        gradI.each_row() -= init.t();
        gradI.each_col() %= init;
        for (unsigned int j = 0; j < emission.n_rows; j++) {
          double tmp = 0.0;
          for (unsigned int r = 0; r < obs.n_rows; r++) {
            tmp += emissionLog(j, obs(r, 0, k), r);
          }
          gradIrow(j) += exp(tmp + beta(j, 0, k) - ll(k));
        }
        
        gradIrow = gradI * gradIrow;
        gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
        countgrad += ind.n_elem;
      }
    }
    grad = sum(gradmat, 1);
  return Rcpp::List::create(Rcpp::Named("objective") = -sum(ll), Rcpp::Named("gradient") = Rcpp::wrap(-grad));
}
