// log-likelihood and gradients of MHMM using log-space
#include "seqHMM.h"
// [[Rcpp::export]]

List log_objectivex(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, IntegerVector transNZ,
    IntegerVector emissNZ, IntegerVector initNZ, const arma::ivec& nSymbols, const arma::mat& coef,
    const arma::mat& X, const arma::ivec& numberOfStates, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, false);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, false);

  arma::imat ANZ(transNZ.begin(), emission.n_rows, emission.n_rows, false);
  arma::icube BNZ(emissNZ.begin(), emission.n_rows, emission.n_cols - 1, emission.n_slices, false);
  arma::ivec INZ(initNZ.begin(), emission.n_rows, false);

  int q = coef.n_rows;
  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem - 1) * q,
      arma::fill::zeros);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    grad.fill(-arma::math::inf());
    return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  }

  weights.each_row() /= sum(weights, 0);
  weights = log(weights);

  arma::vec initLog = log(init);
  arma::mat transitionLog = log(transition);
  arma::cube emissionLog = log(emission);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k



  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = initLog + reparma(weights.col(k), numberOfStates);
  }

  log_internalForwardx(transitionLog, emissionLog, initk, obs, alpha, threads);
  log_internalBackward(transitionLog, emissionLog, obs, beta, threads);

  arma::vec ll(obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
  }

  arma::ivec cumsumstate = arma::cumsum(numberOfStates);

  arma::mat gradmat(
      arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ) + (numberOfStates.n_elem - 1) * q,
      obs.n_slices, arma::fill::zeros);


#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(q, gradmat, nSymbols, ANZ, BNZ, INZ, ll,        \
    numberOfStates, cumsumstate, obs, init, X, weights, transition, emission, \
    initLog, transitionLog, emissionLog, initk, alpha, beta)
  for (int k = 0; k < obs.n_slices; k++) {
    int countgrad = 0;

    // transitionMatrix
    if (arma::accu(ANZ) > 0) {
      for (int jj = 0; jj < numberOfStates.n_elem; jj++) {
        arma::vec gradArow(numberOfStates(jj));
        arma::mat gradA(numberOfStates(jj), numberOfStates(jj));
        int ind_jj = cumsumstate(jj) - numberOfStates(jj);
        for (int i = 0; i < numberOfStates(jj); i++) {
          arma::uvec ind = arma::find(
              ANZ.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1));

          if (ind.n_elem > 0) {
            gradArow.zeros();
            gradA.eye();
            gradA.each_row() -= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1);
            gradA.each_col() %= transition.row(ind_jj + i).subvec(ind_jj, cumsumstate(jj) - 1).t();

            for (int j = 0; j < numberOfStates(jj); j++) {
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                double tmp = 0.0;
                for (unsigned int r = 0; r < obs.n_rows; r++) {
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
                double tmp = 0.0;
                for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                  if (r2 != r) {
                    tmp += emissionLog(i, obs(r2, 0, k), r2);
                  }
                }
                gradBrow(j) += exp(initk(i, k) + tmp + beta(i, 0, k) - ll(k));
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
      for (int i = 0; i < numberOfStates.n_elem; i++) {
        int ind_i = cumsumstate(i) - numberOfStates(i);
        arma::uvec ind = arma::find(INZ.subvec(ind_i, cumsumstate(i) - 1));
        if (ind.n_elem > 0) {
          arma::vec gradIrow(numberOfStates(i), arma::fill::zeros);
          for (int j = 0; j < numberOfStates(i); j++) {
            double tmp = 0.0;
            for (unsigned int r = 0; r < obs.n_rows; r++) {
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
    for (int jj = 1; jj < numberOfStates.n_elem; jj++) {
      int ind_jj = cumsumstate(jj) - numberOfStates(jj);
      for (int j = 0; j < emission.n_rows; j++) {
        double tmp = 0.0;
        for (unsigned int r = 0; r < obs.n_rows; r++) {
          tmp += emissionLog(j, obs(r, 0, k), r);
        }
        if ((j >= ind_jj) & (j < cumsumstate(jj))) {
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
  return List::create(Named("objective") = -sum(ll), Named("gradient") = wrap(-sum(gradmat, 1)));
}
