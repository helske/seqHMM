// log-likelihood and gradients of HMM
#include "pre2-optcoef.h"
#include "pre2-forward_backward.h"
#include "reparma.h"
#include "useomp.h"

// [[Rcpp::export]]
Rcpp::List objective(const arma::mat& transition, const arma::cube& emission,
  const arma::vec& init, arma::ucube& obs, const arma::umat& ANZ,
  const arma::ucube& BNZ, const arma::uvec& INZ, const arma::uvec& nSymbols, arma::uword threads) {
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), arma::fill::zeros);
  arma::uword error = 0;
  double ll = 0;
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) reduction(+:ll) num_threads(threads) \
  default(shared) //shared(grad, nSymbols, ANZ, BNZ, INZ, obs, init, transition, emission, error, arma::fill::zeros)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      if (error == 0) {
        arma::mat alpha(emission.n_rows, obs.n_cols); //m,n
        arma::vec scales(obs.n_cols); //n
        uvForward(transition.t(), emission, init, obs.slice(k), alpha, scales);
        arma::mat beta(emission.n_rows, obs.n_cols); //m,n
        uvBackward(transition, emission, obs.slice(k), beta, scales);
        int countgrad = 0;
        arma::vec grad_k(grad.n_elem, arma::fill::zeros);
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
                double tmp = 1.0;
                for (arma::uword r = 0; r < obs.n_rows; ++r) {
                  tmp *= emission(j, obs(r, t + 1, k), r);
                }
                gradArow(j) += alpha(i, t) * tmp * beta(j, t + 1);
              }

            }
            gradArow = gradA * gradArow;
            grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradArow.rows(ind);
            countgrad += ind.n_elem;
          }
        }
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
                  double tmp = 1.0;
                  for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                    if (r2 != r) {
                      tmp *= emission(i, obs(r2, 0, k), r2);
                    }
                  }
                  gradBrow(j) += init(i) * tmp * beta(i, 0);
                }
                for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                  if (obs(r, t + 1, k) == j) {
                    double tmp = 1.0;
                    for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                      if (r2 != r) {
                        tmp *= emission(i, obs(r2, t + 1, k), r2);
                      }
                    }
                    gradBrow(j) += arma::dot(alpha.col(t), transition.col(i)) * tmp
                      * beta(i, t + 1);
                  }
                }
              }
              gradBrow = gradB * gradBrow;
              grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradBrow.rows(ind);
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
            double tmp = 1.0;
            for (arma::uword r = 0; r < obs.n_rows; ++r) {
              tmp *= emission(j, obs(r, 0, k), r);
            }
            gradIrow(j) += tmp * beta(j, 0);
          }
          gradIrow = gradI * gradIrow;
          grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
          countgrad += ind.n_elem;
        }
        if (!scales.is_finite() || !beta.is_finite()) {
#pragma omp atomic
          error++;
        } else {
          ll -= arma::sum(log(scales));
#pragma omp critical
          grad += grad_k;
        }
      }
    }
    if(error > 0){
      ll = -arma::datum::inf;
      grad.fill(-arma::datum::inf);
    }
    return Rcpp::List::create(Rcpp::Named("objective") = -ll, Rcpp::Named("gradient") = Rcpp::wrap(-grad));
}

// log-likelihood and gradients of MHMM
// [[Rcpp::export]]
Rcpp::List objectivex(const arma::mat& transition, const arma::cube& emission,
                      const arma::vec& init, const arma::ucube& obs, const arma::umat& ANZ,
                      const arma::ucube& BNZ, const arma::uvec& INZ, const arma::uvec& nSymbols,
                      const arma::mat& coef, const arma::mat& X, arma::uvec& numberOfStates,
                      arma::uword threads) {
  
  arma::uword q = coef.n_rows;
  arma::vec grad(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem- 1) * q,
      arma::fill::zeros);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    grad.fill(-arma::datum::inf);
    return Rcpp::List::create(Rcpp::Named("objective") = arma::datum::inf, Rcpp::Named("gradient") = Rcpp::wrap(grad));
  }
  weights.each_row() /= sum(weights, 0);
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }
  arma::uvec cumsumstate = arma::cumsum(numberOfStates);
  arma::uword error = 0;
  double ll = 0;
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) reduction(+:ll) num_threads(threads) default(shared) 
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    if (error == 0) {
      arma::mat alpha(emission.n_rows, obs.n_cols); //m,n
      arma::vec scales(obs.n_cols); //n
      uvForward(transition.t(), emission, initk.col(k), obs.slice(k), alpha, scales);
      arma::mat beta(emission.n_rows, obs.n_cols); //m,n
      uvBackward(transition, emission, obs.slice(k), beta, scales);
      int countgrad = 0;
      arma::vec grad_k(grad.n_elem, arma::fill::zeros);
      // transitionMatrix
      if (arma::accu(ANZ) > 0) {
        for (arma::uword jj = 0; jj < numberOfStates.n_elem; ++jj) {
          arma::vec gradArow(numberOfStates(jj));
          arma::mat gradA(numberOfStates(jj), numberOfStates(jj));
          int ind_jj = cumsumstate(jj) - numberOfStates(jj);
          for (arma::uword i = 0; i < numberOfStates(jj); ++i) {
            arma::uvec ind = arma::find(ANZ.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1));
            if (ind.n_elem > 0) {
              gradArow.zeros();
              gradA.eye();
              gradA.each_row() -= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1);
              gradA.each_col() %= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1).t();
              for (arma::uword j = 0; j < numberOfStates(jj); ++j) {
                for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                  double tmp = alpha(ind_jj + i, t);
                  for (arma::uword r = 0; r < obs.n_rows; ++r) {
                    tmp *= emission(ind_jj + j, obs(r, t + 1, k), r);
                  }
                  gradArow(j) += tmp * beta(ind_jj + j, t + 1);
                }
              }
              gradArow = gradA * gradArow;
              grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradArow.rows(ind);
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
                  double tmp = initk(i, k);
                  for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                    if (r2 != r) {
                      tmp *= emission(i, obs(r2, 0, k), r2);
                    }
                  }
                  gradBrow(j) += tmp * beta(i, 0);
                }
                for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                  if (obs(r, t + 1, k) == j) {
                    double tmp = beta(i, t + 1);
                    for (arma::uword r2 = 0; r2 < obs.n_rows; ++r2) {
                      if (r2 != r) {
                        tmp *= emission(i, obs(r2, t + 1, k), r2);
                      }
                    }
                    gradBrow(j) += arma::dot(alpha.col(t), transition.col(i)) * tmp;
                  }
                }
              }
              gradBrow = gradB * gradBrow;
              grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradBrow.rows(ind);
              countgrad += ind.n_elem;
            }
          }
        }
      }
      if (arma::accu(INZ) > 0) {
        for (arma::uword i = 0; i < numberOfStates.n_elem; ++i) {
          int ind_i = cumsumstate(i) - numberOfStates(i);
          arma::uvec ind = arma::find(
            INZ.subvec(ind_i, cumsumstate(i) - 1));
          if (ind.n_elem > 0) {
            arma::vec gradIrow(numberOfStates(i), arma::fill::zeros);
            for (arma::uword j = 0; j < numberOfStates(i); ++j) {
              double tmp = weights(i, k);
              for (arma::uword r = 0; r < obs.n_rows; ++r) {
                tmp *= emission(ind_i + j, obs(r, 0, k), r);
              }
              gradIrow(j) += tmp * beta(ind_i + j, 0);
            }
            arma::mat gradI(numberOfStates(i), numberOfStates(i), arma::fill::zeros);
            gradI.eye();
            gradI.each_row() -= init.subvec(ind_i, cumsumstate(i) - 1).t();
            gradI.each_col() %= init.subvec(ind_i, cumsumstate(i) - 1);
            gradIrow = gradI * gradIrow;
            grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
            countgrad += ind.n_elem;
          }
        }
      }
      for (arma::uword jj = 1; jj < numberOfStates.n_elem; ++jj) {
        arma::uword ind_jj = (cumsumstate(jj) - numberOfStates(jj));
        
        for (arma::uword j = 0; j < emission.n_rows; ++j) {
          double tmp = 1.0;
          for (arma::uword r = 0; r < obs.n_rows; ++r) {
            tmp *= emission(j, obs(r, 0, k), r);
          }
          if ((j >= ind_jj) && (j < cumsumstate(jj))) {
            grad_k.subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) += tmp
            * beta(j, 0) * initk(j, k) * X.row(k).t() * (1.0 - weights(jj, k));
          } else {
            grad_k.subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) -= tmp
            * beta(j, 0) * initk(j, k) * X.row(k).t() * weights(jj, k);
          }
        }
      }
      if (!scales.is_finite() || !beta.is_finite()) {
#pragma omp atomic
        error++;
      } else {
        ll -= arma::sum(log(scales));
#pragma omp critical
        grad += grad_k;
      }
    }
  }
  if(error > 0){
    ll = -arma::datum::inf;
    grad.fill(-arma::datum::inf);
  }
  return Rcpp::List::create(Rcpp::Named("objective") = -ll, Rcpp::Named("gradient") = Rcpp::wrap(-grad));
}
