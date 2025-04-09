// log-likelihood and gradients of HMM using log-space

#include "pre2-log_forward_backward.h"
#include "logsumexp.h"
#include "reparma.h"
#include "useomp.h"

// [[Rcpp::export]]
Rcpp::List log_objective(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, const arma::ucube& obs, const arma::umat& ANZ,
  const arma::ucube& BNZ, const arma::uvec& INZ, arma::uvec& nSymbols, arma::uword threads) {
  
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), arma::fill::zeros);
  
  arma::vec initLog = log(init);
  arma::mat transitionLog = log(transition);
  arma::cube emissionLog = log(emission);
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  
  log_internalForward(transitionLog, emissionLog, initLog, obs, alpha, threads);
  log_internalBackward(transitionLog, emissionLog, obs, beta, threads);
  
  arma::vec ll(obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
  }
  
  arma::mat gradmat(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), obs.n_slices,
    arma::fill::zeros);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) default(shared)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      int countgrad = 0;
      
      // transitionMatrix
      arma::vec gradArow(emission.n_rows);
      arma::mat gradA(emission.n_rows, emission.n_rows);
      
      for (arma::uword i = 0; i < emission.n_rows; ++i) {
        arma::uvec ind = arma::find(ANZ.row(i));
        
        if (ind.n_elem > 0) {
          gradArow.zeros();
          gradA.eye();
          gradA.each_row() -= transition.row(i);
          gradA.each_col() %= transition.row(i).t();
          
          for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
            for (arma::uword j = 0; j < emission.n_rows; ++j) {
              double tmp = 0.0;
              for (arma::uword r = 0; r < obs.n_rows; ++r) {
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
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        arma::vec gradBrow(nSymbols[r]);
        arma::mat gradB(nSymbols[r], nSymbols[r]);
        for (arma::uword i = 0; i < emission.n_rows; ++i) {
          arma::uvec ind = arma::find(BNZ.slice(r).row(i));
          if (ind.n_elem > 0) {
            gradBrow.zeros();
            gradB.eye();
            gradB.each_row() -= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1);
            gradB.each_col() %= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1).t();
            for (arma::uword j = 0; j < nSymbols[r]; ++j) {
              if (obs(r, 0, k) == j) {
                double tmp = 0.0;
                for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                  if (r2 != r) {
                    tmp += emissionLog(i, obs(r2, 0, k), r2);
                  }
                }
                gradBrow(j) += exp(initLog(i) + tmp + beta(i, 0, k) - ll(k));
              }
              for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                if (obs(r, t + 1, k) == j) {
                  double tmp = 0.0;
                  for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
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
        for (arma::uword j = 0; j < emission.n_rows; ++j) {
          double tmp = 0.0;
          for (arma::uword r = 0; r < obs.n_rows; ++r) {
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

// log-likelihood and gradients of MHMM using log-space
// [[Rcpp::export]]
Rcpp::List log_objectivex(const arma::mat& transition, const arma::cube& emission,
                          const arma::vec& init, const arma::ucube& obs, const arma::umat& ANZ,
                          const arma::ucube& BNZ, const arma::uvec& INZ, const arma::uvec& nSymbols, const arma::mat& coef,
                          const arma::mat& X, const arma::uvec& numberOfStates, arma::uword threads) {
  
  int q = coef.n_rows;
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem - 1) * q,
                 arma::fill::zeros);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    grad.fill(-arma::datum::inf);
    return Rcpp::List::create(Rcpp::Named("objective") = arma::datum::inf, Rcpp::Named("gradient") = Rcpp::wrap(grad));
  }
  
  weights.each_row() /= sum(weights, 0);
  weights = log(weights);
  
  arma::vec initLog = log(init);
  arma::mat transitionLog = log(transition);
  arma::cube emissionLog = log(emission);
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  
  
  
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    initk.col(k) = initLog + reparma(weights.col(k), numberOfStates);
  }
  
  log_internalForward(transitionLog, emissionLog, initk, obs, alpha, threads);
  log_internalBackward(transitionLog, emissionLog, obs, beta, threads);
  
  arma::vec ll(obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
  }
  
  arma::uvec cumsumstate = arma::cumsum(numberOfStates);
  
  arma::mat gradmat(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem - 1) * q,
      obs.n_slices, arma::fill::zeros);
  
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) default(shared)
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    int countgrad = 0;
    
    // transitionMatrix
    if (arma::accu(ANZ) > 0) {
      for (arma::uword jj = 0; jj < numberOfStates.n_elem; ++jj) {
        arma::vec gradArow(numberOfStates(jj));
        arma::mat gradA(numberOfStates(jj), numberOfStates(jj));
        int ind_jj = cumsumstate(jj) - numberOfStates(jj);
        for (arma::uword i = 0; i < numberOfStates(jj); ++i) {
          arma::uvec ind = arma::find(
            ANZ.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1));
          
          if (ind.n_elem > 0) {
            gradArow.zeros();
            gradA.eye();
            gradA.each_row() -= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1);
            gradA.each_col() %= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1).t();
            
            for (arma::uword j = 0; j < numberOfStates(jj); ++j) {
              for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                double tmp = 0.0;
                for (arma::uword r = 0; r < obs.n_rows; ++r) {
                  tmp += emissionLog(ind_jj + j, obs(r, t + 1, k), r);
                }
                gradArow(j) += exp(
                  alpha(ind_jj + i, t, k) + tmp
                + beta(ind_jj + j, t + 1, k) - ll(k));
                
              }
              
            }
            
            gradArow = gradA * gradArow;
            gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradArow.rows(ind);
            countgrad += ind.n_elem;
          }
        }
      }
    }
    if (arma::accu(BNZ) > 0) {
      // emissionMatrix
      for (arma::uword r = 0; r < obs.n_rows; ++r) {
        arma::vec gradBrow(nSymbols(r));
        arma::mat gradB(nSymbols(r), nSymbols(r));
        for (arma::uword i = 0; i < emission.n_rows; ++i) {
          arma::uvec ind = arma::find(BNZ.slice(r).row(i));
          if (ind.n_elem > 0) {
            gradBrow.zeros();
            gradB.eye();
            gradB.each_row() -= emission.slice(r).row(i).subvec(0, nSymbols(r) - 1);
            gradB.each_col() %= emission.slice(r).row(i).subvec(0, nSymbols(r) - 1).t();
            for (arma::uword j = 0; j < nSymbols(r); ++j) {
              if (obs(r, 0, k) == j) {
                double tmp = 0.0;
                for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                  if (r2 != r) {
                    tmp += emissionLog(i, obs(r2, 0, k), r2);
                  }
                }
                gradBrow(j) += exp(initk(i, k) + tmp + beta(i, 0, k) - ll(k));
              }
              for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                if (obs(r, t + 1, k) == j) {
                  double tmp = 0.0;
                  for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                    if (r2 != r) {
                      tmp += emissionLog(i, obs(r2, t + 1, k), r2);
                    }
                  }
                  gradBrow(j) += arma::accu(
                    exp(alpha.slice(k).col(t) + tmp + transitionLog.col(i) + beta(i, t + 1, k) - ll(k)));
                }
              }
              
            }
            gradBrow = gradB * gradBrow;
            gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradBrow.rows(ind);
            countgrad += ind.n_elem;
            
          }
        }
      }
    }
    if (arma::accu(INZ) > 0) {
      for (arma::uword i = 0; i < numberOfStates.n_elem; ++i) {
        int ind_i = cumsumstate(i) - numberOfStates(i);
        arma::uvec ind = arma::find(INZ.subvec(ind_i, cumsumstate(i) - 1));
        if (ind.n_elem > 0) {
          arma::vec gradIrow(numberOfStates(i), arma::fill::zeros);
          for (arma::uword j = 0; j < numberOfStates(i); ++j) {
            double tmp = 0.0;
            for (arma::uword r = 0; r < obs.n_rows; ++r) {
              tmp += emissionLog(ind_i + j, obs(r, 0, k), r);
            }
            gradIrow(j) += exp(tmp + beta(ind_i + j, 0, k) - ll(k) + weights(i, k));
            
          }
          arma::mat gradI(numberOfStates(i), numberOfStates(i), arma::fill::zeros);
          gradI.eye();
          gradI.each_row() -= init.subvec(ind_i, cumsumstate(i) - 1).t();
          gradI.each_col() %= init.subvec(ind_i, cumsumstate(i) - 1);
          gradIrow = gradI * gradIrow;
          gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
          countgrad += ind.n_elem;
        }
      }
    }
    for (arma::uword jj = 1; jj < numberOfStates.n_elem; ++jj) {
      arma::uword ind_jj = cumsumstate(jj) - numberOfStates(jj);
      for (arma::uword j = 0; j < emission.n_rows; ++j) {
        double tmp = 0.0;
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          tmp += emissionLog(j, obs(r, 0, k), r);
        }
        if ((j >= ind_jj) && (j < cumsumstate(jj))) {
          gradmat.col(k).subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) += exp(
            tmp + beta(j, 0, k) - ll(k) + initk(j, k)) * X.row(k).t()
          * (1.0 - exp(weights(jj, k)));
        } else {
          gradmat.col(k).subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) -= exp(
            tmp + beta(j, 0, k) - ll(k) + initk(j, k)) * X.row(k).t() * exp(weights(jj, k));
        }
      }
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("objective") = -sum(ll), Rcpp::Named("gradient") = Rcpp::wrap(-sum(gradmat, 1)));
}
