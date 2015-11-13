#include "seqHMM.h"
// [[Rcpp::export]]

List log_EMx(NumericVector transitionMatrix, NumericVector emissionArray,
    NumericVector initialProbs, IntegerVector obsArray, IntegerVector nSymbols, NumericMatrix coefs,
    NumericMatrix X_, IntegerVector numberOfStates, int itermax, double tol, int trace, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);

  int q = coefs.nrow();
  arma::mat coef(coefs.begin(), q, coefs.ncol());
  coef.col(0).zeros();
  arma::mat X(X_.begin(), obs.n_rows, q, false);
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return List::create(Named("error") = 1);
  }
  weights.each_row() /= sum(weights, 0);
  weights = log(weights);

  transition = log(transition);
  emission = log(emission);
  init = log(init);

  arma::mat initk(emission.n_rows, obs.n_rows);
  for (unsigned int k = 0; k < obs.n_rows; k++) {
    initk.col(k) = init + reparma(weights.col(k), numberOfStates);
  }

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_rows); //m,n,k

  log_internalForwardx(transition, emission, initk, obs, alpha, threads);
  log_internalBackward(transition, emission, obs, beta, threads);

  arma::vec ll(obs.n_rows);

  for (int k = 0; k < oDims[0]; k++) {
    ll(k) = logSumExp(alpha.slice(k).col(oDims[1] - 1));
  }

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

    arma::mat ksii(emission.n_rows, emission.n_rows, arma::fill::zeros);
    arma::cube gamma(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
    arma::vec delta(emission.n_rows, arma::fill::zeros);

    for (unsigned int k = 0; k < obs.n_rows; k++) {
      delta += exp(alpha.slice(k).col(0) + beta.slice(k).col(0) - ll(k));
    }

#pragma omp parallel for if(obs.n_rows>=threads) schedule(static) num_threads(threads) \
    default(none) shared(transition, obs, ll, alpha, beta, emission, ksii, gamma, nSymbols)
    for (int k = 0; k < obs.n_rows; k++) {
      if (obs.n_cols > 1) {
        for (unsigned int j = 0; j < emission.n_rows; j++) {
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            if (transition(i, j) > -arma::math::inf()) {
              arma::vec tmpnm1(obs.n_cols - 1);
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                tmpnm1(t) = alpha(i, t, k) + transition(i, j) + beta(j, t + 1, k);
                for (unsigned int r = 0; r < obs.n_slices; r++) {
                  tmpnm1(t) += emission(j, obs(k, t + 1, r), r);
                }
              }
#pragma omp atomic
              ksii(i, j) += exp(logSumExp(tmpnm1) - ll(k));
            }
          }
        }
      }

      for (unsigned int r = 0; r < emission.n_slices; r++) {
        for (int l = 0; l < nSymbols[r]; l++) {
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            if (emission(i, l, r) > -arma::math::inf()) {
              arma::vec tmpn(obs.n_cols);
              for (unsigned int t = 0; t < obs.n_cols; t++) {
                if (l == (obs(k, t, r))) {
                  tmpn(t) = alpha(i, t, k) + beta(i, t, k);
                } else
                  tmpn(t) = -arma::math::inf();
              }
#pragma omp atomic
              gamma(i, l, r) += exp(logSumExp(tmpn) - ll(k));

            }
          }
        }
      }
    }

    unsigned int error = log_optCoef(weights, obs, emission, initk, beta, ll, coef, X, cumsumstate,
        numberOfStates, trace);
    if (error != 0) {
      return List::create(Named("error") = error);
    }
    if (obs.n_cols > 1) {
      ksii.each_col() /= sum(ksii, 1);
      transition = log(ksii);
    }
    for (unsigned int r = 0; r < emission.n_slices; r++) {
      gamma.slice(r).cols(0, nSymbols(r) - 1).each_col() /= sum(
          gamma.slice(r).cols(0, nSymbols(r) - 1), 1);
      emission.slice(r).cols(0, nSymbols(r) - 1) = log(gamma.slice(r).cols(0, nSymbols(r) - 1));
    }

    for (int i = 0; i < numberOfStates.size(); i++) {
      delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1) /= arma::as_scalar(
          arma::accu(delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1)));
    }

    init = log(delta);

    for (unsigned int k = 0; k < obs.n_rows; k++) {
      initk.col(k) = init + reparma(weights.col(k), numberOfStates);
    }

    log_internalForwardx(transition, emission, initk, obs, alpha, threads);
    log_internalBackward(transition, emission, obs, beta, threads);

    for (int k = 0; k < oDims[0]; k++) {
      ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
    }

    double tmp = sum(ll);
    change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
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
  return List::create(Named("coefficients") = wrap(coef), Named("initialProbs") = wrap(exp(init)),
      Named("transitionMatrix") = wrap(exp(transition)),
      Named("emissionArray") = wrap(exp(emission)), Named("logLik") = sumlogLik,
      Named("iterations") = iter, Named("change") = change, Named("error") = 0);
}
