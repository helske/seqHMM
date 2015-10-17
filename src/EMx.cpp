#include "seqHMM.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List EMx(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
    IntegerVector obsArray, IntegerVector nSymbols, NumericMatrix coefs, NumericMatrix X_,
    IntegerVector numberOfStates, int itermax = 100, double tol = 1e-8, int trace = 0, int threads =
        1) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::vec init(initialProbs.begin(), eDims[0], true);
  arma::mat transition(transitionMatrix.begin(), eDims[0], eDims[0], true);
  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);

  int q = coefs.nrow();
  arma::mat coef(coefs.begin(), q, coefs.ncol());
  coef.col(0).zeros();
  arma::mat X(X_.begin(), oDims[0], q);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return List::create(Named("error") = 1);
  }
  weights.each_row() /= sum(weights, 0);

  arma::mat initk(eDims[0], oDims[0]);
  for (int k = 0; k < oDims[0]; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }

  arma::cube alpha(eDims[0], oDims[1], oDims[0]); //m,n,k
  arma::cube beta(eDims[0], oDims[1], oDims[0]); //m,n,k
  arma::mat scales(oDims[1], oDims[0]); //m,n,k

  internalForwardx(transition, emission, initk, obs, alpha, scales, threads);
  internalBackward(transition, emission, obs, beta, scales, threads);

  arma::rowvec ll = arma::sum(log(scales));
  double sumlogLik = sum(ll);

  if (trace > 0) {
    Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
  }
  //  
  //  //EM-algorithm begins
  //  
  double change = tol + 1.0;
  int iter = 0;

  IntegerVector cumsumstate = cumsum(numberOfStates);

  while ((change > tol) & (iter < itermax)) {
    iter++;

    arma::mat ksii(eDims[0], eDims[0], arma::fill::zeros);
    arma::cube gamma(eDims[0], eDims[1], eDims[2], arma::fill::zeros);
    arma::vec delta(eDims[0], arma::fill::zeros);

    for (int k = 0; k < oDims[0]; k++) {
      delta += alpha.slice(k).col(0) % beta.slice(k).col(0);
    }
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int k = 0; k < oDims[0]; k++) {
      for (int i = 0; i < eDims[0]; i++) {
        for (int j = 0; j < eDims[0]; j++) {
          if (transition(i, j) > 0.0) {
            for (int t = 0; t < (oDims[1] - 1); t++) {
              double tmp = alpha(i, t, k) * transition(i, j) * beta(j, t + 1, k) / scales(t + 1, k);
              for (int r = 0; r < oDims[2]; r++) {
                tmp *= emission(j, obs(k, t + 1, r), r);
              }
#pragma omp atomic
              ksii(i, j) += tmp;
            }
          }
        }
      }

      for (int r = 0; r < eDims[2]; r++) {
        for (int i = 0; i < eDims[0]; i++) {
          for (int l = 0; l < nSymbols[r]; l++) {
            if (emission(i, l, r) > 0.0) {
              for (int t = 0; t < oDims[1]; t++) {
                if (l == (obs(k, t, r))) {
                  double tmp = alpha(i, t, k) * beta(i, t, k);
#pragma omp atomic
                  gamma(i, l, r) += tmp;
                }
              }
            }
          }
        }
      }
    }

    unsigned int error = optCoef(weights, obs, emission, initk, beta, scales, coef, X, cumsumstate,
        numberOfStates, trace);
    if (error != 0) {
      return List::create(Named("error") = error);
    }

    if (oDims[1] > 1) {
      ksii.each_col() /= sum(ksii, 1);
      transition = ksii;
    }
    for (int r = 0; r < eDims[2]; r++) {

      gamma.slice(r).cols(0, nSymbols(r) - 1).each_col() /= sum(
          gamma.slice(r).cols(0, nSymbols(r) - 1), 1);
      emission.slice(r).cols(0, nSymbols(r) - 1) = gamma.slice(r).cols(0, nSymbols(r) - 1);
    }

    for (int i = 0; i < numberOfStates.size(); i++) {
      delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1) /= arma::as_scalar(
          arma::accu(delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1)));
    }
    init = delta;

    for (int k = 0; k < oDims[0]; k++) {
      initk.col(k) = init % reparma(weights.col(k), numberOfStates);
    }

    internalForwardx(transition, emission, initk, obs, alpha, scales, threads);
    internalBackward(transition, emission, obs, beta, scales, threads);

    ll = sum(log(scales));
    double tmp = sum(ll);
    change = (tmp - sumlogLik) / (abs(sumlogLik) + 0.1);
    sumlogLik = tmp;
    if (!arma::is_finite(sumlogLik)) {
      return List::create(Named("error") = 4);
    }
    if (trace > 1) {
      Rcout << "iter: " << iter;
      Rcout << " logLik: " << sumlogLik;
      Rcout << " relative change: " << change << std::endl;
    }

  }
  if (trace > 0) {
    if (iter == itermax) {
      Rcpp::Rcout << "EM algorithm stopped after reaching the maximum number of " << iter
          << " iterations." << std::endl;
    } else {
      Rcpp::Rcout << "EM algorithm stopped after reaching the relative change of " << change;
      Rcpp::Rcout << " after " << iter << " iterations." << std::endl;
    }
    Rcpp::Rcout << "Final log-likelihood: " << sumlogLik << std::endl;
  }
  return List::create(Named("coefficients") = wrap(coef), Named("initialProbs") = wrap(init),
      Named("transitionMatrix") = wrap(transition), Named("emissionArray") = wrap(emission),
      Named("logLik") = sumlogLik, Named("iterations") = iter, Named("change") = change,
      Named("error") = 0);
}
