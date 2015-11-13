#include "seqHMM.h"
// [[Rcpp::export]]

List objectivex(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, IntegerVector transNZ,
    IntegerVector emissNZ, IntegerVector initNZ, IntegerVector nSymbols, NumericMatrix coefs,
    NumericMatrix X_, IntegerVector numberOfStates, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);

  arma::imat ANZ(transNZ.begin(), emission.n_rows, emission.n_rows, false);
  arma::icube BNZ(emissNZ.begin(), emission.n_rows, emission.n_cols - 1, emission.n_slices, false);
  arma::ivec INZ(initNZ.begin(), emission.n_rows, false);

  int q = coefs.nrow();
  arma::vec grad(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.size() - 1) * q,
      arma::fill::zeros);
  arma::mat coef(coefs.begin(), q, numberOfStates.size());
  coef.col(0).zeros();
  arma::mat X(X_.begin(), obs.n_rows, q);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad),
        Named("error") = 1);
  }

  weights.each_row() /= sum(weights, 0);

  arma::mat initk(emission.n_rows, obs.n_rows);

  for (unsigned int k = 0; k < obs.n_rows; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k
  arma::mat scales(obs.n_cols, obs.n_rows); //m,n,k

  internalForwardx(transition, emission, initk, obs, alpha, scales, threads);
  if (!alpha.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  }
  internalBackward(transition, emission, obs, beta, scales, threads);
  if (!beta.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad),
        Named("error") = 2);
  }

  IntegerVector cumsumstate = cumsum(numberOfStates);

  arma::mat gradmat(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.size() - 1) * q,
      obs.n_rows, arma::fill::zeros);

#pragma omp parallel for if(obs.n_rows >= threads) schedule(static) num_threads(threads) \
default(none) shared(q, alpha, beta, scales, gradmat, nSymbols, ANZ, BNZ, INZ, \
  numberOfStates, cumsumstate, obs, init, initk, X, weights, transition, emission)
  for (int k = 0; k < obs.n_rows; k++) {
    int countgrad = 0;
    // transitionMatrix
    if (arma::accu(ANZ) > 0) {
      for (int jj = 0; jj < numberOfStates.size(); jj++) {
        arma::vec gradArow(numberOfStates(jj));
        arma::mat gradA(numberOfStates(jj), numberOfStates(jj));
        for (int i = 0; i < numberOfStates(jj); i++) {
          arma::uvec ind = arma::find(
              ANZ.row(cumsumstate(jj) - numberOfStates(jj) + i).subvec(
                  cumsumstate(jj) - numberOfStates(jj), cumsumstate(jj) - 1));

          if (ind.n_elem > 0) {
            gradArow.zeros();
            gradA.eye();
            gradA.each_row() -= transition.row(cumsumstate(jj) - numberOfStates(jj) + i).subvec(
                cumsumstate(jj) - numberOfStates(jj), cumsumstate(jj) - 1);
            gradA.each_col() %= transition.row(cumsumstate(jj) - numberOfStates(jj) + i).subvec(
                cumsumstate(jj) - numberOfStates(jj), cumsumstate(jj) - 1).t();

            for (int j = 0; j < numberOfStates(jj); j++) {
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                double tmp = 1.0;
                for (unsigned int r = 0; r < obs.n_slices; r++) {
                  tmp *= emission(cumsumstate(jj) - numberOfStates(jj) + j, obs(k, t + 1, r), r);
                }
                gradArow(j) += alpha(cumsumstate(jj) - numberOfStates(jj) + i, t, k) * tmp
                    * beta(cumsumstate(jj) - numberOfStates(jj) + j, t + 1, k) / scales(t + 1, k);
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
      for (unsigned int r = 0; r < obs.n_slices; r++) {
        arma::vec gradBrow(nSymbols[r]);
        arma::mat gradB(nSymbols[r], nSymbols[r]);
        for (unsigned int i = 0; i < emission.n_rows; i++) {
          arma::uvec ind = arma::find(BNZ.slice(r).row(i));
          if (ind.n_elem > 0) {
            gradBrow.zeros();
            gradB.eye();
            gradB.each_row() -= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1);
            gradB.each_col() %= emission.slice(r).row(i).subvec(0, nSymbols[r] - 1).t();
            for (int j = 0; j < nSymbols[r]; j++) {
              if (obs(k, 0, r) == j) {
                double tmp = 1.0;
                for (unsigned int r2 = 0; r2 < obs.n_slices; r2++) {
                  if (r2 != r) {
                    tmp *= emission(i, obs(k, 0, r2), r2);
                  }
                }
                gradBrow(j) += initk(i, k) * tmp * beta(i, 0, k) / scales(0, k);
              }
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                if (obs(k, t + 1, r) == j) {
                  double tmp = 1.0;
                  for (unsigned int r2 = 0; r2 < obs.n_slices; r2++) {
                    if (r2 != r) {
                      tmp *= emission(i, obs(k, t + 1, r2), r2);
                    }
                  }
                  gradBrow(j) += arma::dot(alpha.slice(k).col(t), transition.col(i)) * tmp
                      * beta(i, t + 1, k) / scales(t + 1, k);
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
      for (int i = 0; i < numberOfStates.size(); i++) {
        arma::uvec ind = arma::find(
            INZ.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1));
        if (ind.n_elem > 0) {
          arma::vec gradIrow(numberOfStates(i), arma::fill::zeros);
          for (int j = 0; j < numberOfStates(i); j++) {
            double tmp = 1.0;
            for (unsigned int r = 0; r < obs.n_slices; r++) {
              tmp *= emission(cumsumstate(i) - numberOfStates(i) + j, obs(k, 0, r), r);
            }
            gradIrow(j) += tmp * beta(cumsumstate(i) - numberOfStates(i) + j, 0, k) / scales(0, k)
                * weights(i, k);

          }
          arma::mat gradI(numberOfStates(i), numberOfStates(i), arma::fill::zeros);
          gradI.eye();
          gradI.each_row() -=
              init.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1).t();
          gradI.each_col() %= init.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1);
          gradIrow = gradI * gradIrow;
          gradmat.col(k).subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
          countgrad += ind.n_elem;
        }
      }
    }
    for (int jj = 1; jj < numberOfStates.size(); jj++) {
      for (int j = 0; j < emission.n_rows; j++) {
        double tmp = 1.0;
        for (unsigned int r = 0; r < obs.n_slices; r++) {
          tmp *= emission(j, obs(k, 0, r), r);
        }
        if ((j >= (cumsumstate(jj) - numberOfStates(jj))) & (j < cumsumstate(jj))) {
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
      Named("gradient") = wrap(-sum(gradmat, 1)), Named("error") = 0);
}
