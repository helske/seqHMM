// log-likelihood and gradients of HMM
#include "seqHMM.h"

// [[Rcpp::export]]

List objective(const arma::mat& transition, NumericVector emissionArray,
  const arma::vec& init, IntegerVector obsArray, const arma::imat& ANZ,
  IntegerVector emissNZ, const arma::ivec& INZ, const arma::ivec& nSymbols, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], false, true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  arma::icube BNZ(emissNZ.begin(), emission.n_rows, emission.n_cols - 1, emission.n_slices, false, true);

  arma::vec grad(arma::accu(ANZ) + arma::accu(BNZ) + arma::accu(INZ), arma::fill::zeros);

  // arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  // arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  // arma::mat scales(obs.n_cols, obs.n_slices); //m,n,k
  //
  // internalForward(transition, emission, init, obs, alpha, scales, threads);
  // if (!scales.is_finite()) {
  //   grad.fill(-arma::math::inf());
  //   return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  // }
  //
  // internalBackward(transition, emission, obs, beta, scales, threads);
  // if (!beta.is_finite()) {
  //   grad.fill(-arma::math::inf());
  //   return List::create(Named("objective") = arma::math::inf(), Named("gradient") = wrap(grad));
  // }


  unsigned int error = 0;
  double ll = 0;
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) reduction(+:ll) num_threads(threads) \
  default(none) shared(grad, nSymbols, ANZ, BNZ, INZ,                                                      \
    obs, init, transition, emission, error)
    for (int k = 0; k < obs.n_slices; k++) {
      if (error == 0) {
        arma::mat alpha(emission.n_rows, obs.n_cols); //m,n
        arma::mat beta(emission.n_rows, obs.n_cols); //m,n
        arma::vec scales(obs.n_cols); //n

        uvForward(transition, emission, init, obs.slice(k), alpha, scales);
        uvBackward(transition, emission, obs.slice(k), beta, scales);

        int countgrad = 0;
        arma::vec grad_k(grad.n_elem, arma::fill::zeros);
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
                double tmp = 1.0;
                for (unsigned int r = 0; r < obs.n_rows; r++) {
                  tmp *= emission(j, obs(r, t + 1, k), r);
                }
                gradArow(j) += alpha(i, t) * tmp * beta(j, t + 1) / scales(t + 1);
              }

            }

            gradArow = gradA * gradArow;
            grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradArow.rows(ind);
            countgrad += ind.n_elem;
          }
        }
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
                  double tmp = 1.0;
                  for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                    if (r2 != r) {
                      tmp *= emission(i, obs(r2, 0, k), r2);
                    }
                  }
                  gradBrow(j) += init(i) * tmp * beta(i, 0) / scales(0);
                }
                for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                  if (obs(r, t + 1, k) == j) {
                    double tmp = 1.0;
                    for (unsigned int r2 = 0; r2 < obs.n_rows; r2++) {
                      if (r2 != r) {
                        tmp *= emission(i, obs(r2, t + 1, k), r2);
                      }
                    }
                    gradBrow(j) += arma::dot(alpha.col(t), transition.col(i)) * tmp
                      * beta(i, t + 1) / scales(t + 1);
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
          for (unsigned int j = 0; j < emission.n_rows; j++) {
            double tmp = 1.0;
            for (unsigned int r = 0; r < obs.n_rows; r++) {
              tmp *= emission(j, obs(r, 0, k), r);
            }
            gradIrow(j) += tmp * beta(j, 0) / scales(0);
          }

          gradIrow = gradI * gradIrow;
          grad_k.subvec(countgrad, countgrad + ind.n_elem - 1) = gradIrow.rows(ind);
          countgrad += ind.n_elem;
        }
        if (!scales.is_finite() || !beta.is_finite()) {
#pragma omp atomic
          error++;
        } else {
          ll += arma::sum(log(scales));
#pragma omp critical
          grad += grad_k;
//           for (unsigned int ii = 0; ii < grad_k.n_elem; ii++) {
// #pragma omp atomic
//             grad(ii) += grad_k(ii);
//         }
        }
      }
    }
    if(error > 0){
      ll = -arma::math::inf();
      grad.fill(-arma::math::inf());
    }
    return List::create(Named("objective") = -ll, Named("gradient") = wrap(-grad));
}
