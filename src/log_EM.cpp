// EM algorithm for non-mixture hidden Markov models using logarithm space

#include "seqHMM.h"
// [[Rcpp::export]]

List log_EM(NumericVector transitionMatrix, NumericVector emissionArray, NumericVector initialProbs,
    IntegerVector obsArray, const arma::ivec& nSymbols, int itermax, double tol, int trace, int threads) {

  IntegerVector eDims = emissionArray.attr("dim"); //m,p,r
  IntegerVector oDims = obsArray.attr("dim"); //k,n,r

  arma::cube emission(emissionArray.begin(), eDims[0], eDims[1], eDims[2], true);
  arma::icube obs(obsArray.begin(), oDims[0], oDims[1], oDims[2], false, true);
  arma::vec init(initialProbs.begin(), emission.n_rows, true);
  arma::mat transition(transitionMatrix.begin(), emission.n_rows, emission.n_rows, true);

  transition = log(transition);
  emission = log(emission);
  init = log(init);

  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k

  log_internalForward(transition, emission, init, obs, alpha, threads);
  log_internalBackward(transition, emission, obs, beta, threads);

  arma::vec ll(obs.n_slices);

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(obs, alpha, ll)
  for (int k = 0; k < obs.n_slices; k++) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
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

  while ((change > tol) & (iter < itermax)) {
    iter++;

    arma::mat ksii(emission.n_rows, emission.n_rows, arma::fill::zeros);
    arma::cube gamma(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
    arma::vec delta(emission.n_rows, arma::fill::zeros);

    for (unsigned int k = 0; k < obs.n_slices; k++) {
      delta += exp(alpha.slice(k).col(0) + beta.slice(k).col(0) - ll(k));
    }

#pragma omp parallel for if(obs.n_slices>=threads) schedule(static) num_threads(threads) \
    default(none) shared(transition, obs, alpha, beta, ll,                  \
      emission, ksii, gamma, nSymbols)
    for (int k = 0; k < obs.n_slices; k++) {
      if (obs.n_cols > 1) {
        for (unsigned int j = 0; j < emission.n_rows; j++) {
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            if (transition(i, j) > -arma::math::inf()) {
              arma::vec tmpnm1(obs.n_cols - 1);
              for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                tmpnm1(t) = alpha(i, t, k) + transition(i, j) + beta(j, t + 1, k);
                for (unsigned int r = 0; r < obs.n_rows; r++) {
                  tmpnm1(t) += emission(j, obs(r, t + 1, k), r);
                }
              }
#pragma omp atomic
              ksii(i, j) += exp(logSumExp(tmpnm1) - ll(k));
            }
          }
        }
      }

      for (unsigned int r = 0; r < emission.n_slices; r++) {
        for (int l = 0; l < nSymbols(r); l++) {
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            if (emission(i, l, r) > -arma::math::inf()) {
              arma::vec tmpn(obs.n_cols);
              for (unsigned int t = 0; t < obs.n_cols; t++) {
                if (l == (obs(r, t, k))) {
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

    if (obs.n_cols > 1) {
      ksii.each_col() /= sum(ksii, 1);
      transition = log(ksii);
    }
    for (unsigned int r = 0; r < emission.n_slices; r++) {
      gamma.slice(r).cols(0, nSymbols(r) - 1).each_col() /= sum(
          gamma.slice(r).cols(0, nSymbols(r) - 1), 1);
      emission.slice(r).cols(0, nSymbols(r) - 1) = log(gamma.slice(r).cols(0, nSymbols(r) - 1));
    }

    delta /= arma::as_scalar(arma::accu(delta));

    init = log(delta);

    log_internalForward(transition, emission, init, obs, alpha, threads);
    log_internalBackward(transition, emission, obs, beta, threads);
    
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
    default(none) shared(obs, alpha, ll)
    for (int k = 0; k < obs.n_slices; k++) {
      ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
    }

    double tmp = sum(ll);
    change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
    sumlogLik = tmp;
    if (!arma::is_finite(sumlogLik)) {
      return List::create(Named("error") = 6);
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
  return List::create(Named("initialProbs") = wrap(exp(init)),
      Named("transitionMatrix") = wrap(exp(transition)),
      Named("emissionArray") = wrap(exp(emission)), Named("logLik") = sumlogLik,
      Named("iterations") = iter, Named("change") = change);
}
