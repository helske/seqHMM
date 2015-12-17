// log-likelihood and gradients of MHMM
#include "seqHMM.h"
// [[Rcpp::export]]

List objectivex(const arma::mat& transition, NumericVector emissionArray,
                const arma::vec& init, IntegerVector obsArray, const arma::imat& ANZ,
                IntegerVector emissNZ, const arma::ivec& INZ, const arma::ivec& nSymbols,
                const arma::mat& coef, const arma::mat& X, arma::ivec& numberOfStates,
                int threads) {


  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);

  arma::icube BNZ(emissNZ.begin(), emission.n_rows, emission.n_cols - 1, emission.n_slices, false);

  unsigned int q = coef.n_rows;
  arma::vec grad(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem- 1) * q,
      arma::fill::zeros);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  }

  weights.each_row() /= sum(weights, 0);

  arma::mat initk(emission.n_rows, obs.n_slices);

  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_slices); //m,n,k

  arma::sp_mat sp_trans(transition);
  internalForwardx(sp_trans.t(), emission, initk, obs, alpha, scales, threads);
  if (!scales.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  }

  internalBackwardx(sp_trans, emission, obs, beta, scales, threads);
  if (!beta.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  }

  arma::ivec cumsumstate = arma::cumsum(numberOfStates);

  arma::mat gradmat(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem- 1) * q,
      obs.n_slices, arma::fill::zeros);

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads)       \
  default(none) shared(q, alpha, beta, scales, gradmat, nSymbols, ANZ, BNZ, INZ,          \
          numberOfStates, cumsumstate, obs, init, initk, X, weights, transition, emission)
    for (int k = 0; k < obs.n_slices; k++) {
      int countgrad = 0;
      // transitionMatrix
      if (arma::accu(ANZ) > 0) {

        for (int jj = 0; jj < numberOfStates.n_elem; jj++) {
          arma::vec gradArow(numberOfStates(jj));
          arma::mat gradA(numberOfStates(jj), numberOfStates(jj));
          int ind_jj = cumsumstate(jj) - numberOfStates(jj);

          for (int i = 0; i < numberOfStates(jj); i++) {
            arma::uvec ind = arma::find(ANZ.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1));

            if (ind.n_elem > 0) {
              gradArow.zeros();
              gradA.eye();
              gradA.each_row() -= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1);
              gradA.each_col() %= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1).t();


              for (int j = 0; j < numberOfStates(jj); j++) {
                for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                  double tmp = alpha(ind_jj + i, t, k);
                  for (unsigned int r = 0; r < obs.n_rows; r++) {
                    tmp *= emission(ind_jj + j, obs(r, t + 1, k), r);
                  }
                  gradArow(j) += tmp * beta(ind_jj + j, t + 1, k) / scales(t + 1, k);
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
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          arma::vec gradBrow(nSymbols(r));
          arma::mat gradB(nSymbols(r), nSymbols(r));
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            arma::uvec ind = arma::find(BNZ.slice(r).row(i));
            if (ind.n_elem > 0) {
              gradBrow.zeros();
              gradB.eye();
              gradB.each_row() -= emission.slice(r).row(i).subvec(0, nSymbols(r) - 1);
              gradB.each_col() %= emission.slice(r).row(i).subvec(0, nSymbols(r) - 1).t();
              for (int j = 0; j < nSymbols(r); j++) {
                if (obs(r, 0, k) == j) {
                  double tmp = initk(i, k);
                  for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                    if (r2 != r) {
                      tmp *= emission(i, obs(r2, 0, k), r2);
                    }
                  }
                  gradBrow(j) += tmp * beta(i, 0, k) / scales(0, k);
                }
                for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                  if (obs(r, t + 1, k) == j) {
                    double tmp = beta(i, t + 1, k) / scales(t + 1, k);
                    for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                      if (r2 != r) {
                        tmp *= emission(i, obs(r2, t + 1, k), r2);
                      }
                    }
                    gradBrow(j) += arma::dot(alpha.slice(k).col(t), transition.col(i)) * tmp;
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
        for (int i = 0; i < numberOfStates.n_elem; i++) {
          int ind_i = cumsumstate(i) - numberOfStates(i);
          arma::uvec ind = arma::find(
            INZ.subvec(ind_i, cumsumstate(i) - 1));
          if (ind.n_elem > 0) {
            arma::vec gradIrow(numberOfStates(i), arma::fill::zeros);
            for (int j = 0; j < numberOfStates(i); j++) {
              double tmp = weights(i, k);
              for (unsigned int r = 0; r < obs.n_rows; r++) {
                tmp *= emission(ind_i + j, obs(r, 0, k), r);
              }
              gradIrow(j) += tmp * beta(ind_i + j, 0, k) / scales(0, k);

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
      for (int jj = 1; jj < numberOfStates.n_elem; jj++) {
        int ind_jj = (cumsumstate(jj) - numberOfStates(jj));

        for (int j = 0; j < emission.n_rows; j++) {
          double tmp = 1.0;
          for (unsigned int r = 0; r < obs.n_rows; r++) {
            tmp *= emission(j, obs(r, 0, k), r);
          }
          if ((j >= ind_jj) & (j < cumsumstate(jj))) {
            gradmat.col(k).subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) += tmp
            * beta(j, 0, k) / scales(0, k) * initk(j, k) * X.row(k).t() * (1.0 - weights(jj, k));
          } else {
            gradmat.col(k).subvec(countgrad + q * (jj - 1), countgrad + q * jj - 1) -= tmp
            * beta(j, 0, k) / scales(0, k) * initk(j, k) * X.row(k).t() * weights(jj, k);
          }
        }

      }
    }
    return List::create(Named("objective") = -arma::accu(log(scales)),
                        Named("gradient") = wrap(-sum(gradmat, 1)));
}
