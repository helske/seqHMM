// EM algorithm for non-mixture hidden Markov models using logarithm space

#include "pre2-log_forward_backward.h"
#include "pre2-optcoef.h"
#include "logsumexp.h"
#include "reparma.h"
#include "useomp.h"

// [[Rcpp::export]]
Rcpp::List log_EM(const arma::mat& transition_, const arma::cube& emission_, 
  const arma::vec& init_, const arma::ucube& obs, const arma::uvec& nSymbols, 
  int itermax, double tol, int trace, arma::uword threads) {

  // Make sure we don't alter the original vec/mat/cube
  // needed for cube, in future maybe in other cases as well
  arma::cube emission = log(emission_);
  arma::mat transition = log(transition_);
  arma::vec init = log(init_);
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k

  log_internalForward(transition, emission, init, obs, alpha, threads);
  log_internalBackward(transition, emission, obs, beta, threads);

  arma::vec ll(obs.n_slices);

#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(obs, alpha, ll)
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
  }

  double sumlogLik = sum(ll);
  if (trace > 0) {
    Rcpp::Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
  }
  //  
  //  //EM-algorithm begins
  //  
  double change = tol + 1.0;
  int iter = 0;

  while ((change > tol) && (iter < itermax)) {
    iter++;

    arma::mat ksii(emission.n_rows, emission.n_rows, arma::fill::zeros);
    arma::cube gamma(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
    arma::vec delta(emission.n_rows, arma::fill::zeros);

    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      delta += exp(alpha.slice(k).col(0) + beta.slice(k).col(0) - ll(k));
    }

#pragma omp parallel for if(obs.n_slices>=threads) schedule(static) num_threads(threads) \
    default(none) shared(transition, obs, alpha, beta, ll,                  \
      emission, ksii, gamma, nSymbols)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      if (obs.n_cols > 1) {
        for (arma::uword j = 0; j < emission.n_rows; ++j) {
          for (arma::uword i = 0; i < emission.n_rows; ++i) {
            if (transition(i, j) > -arma::datum::inf) {
              arma::vec tmpnm1(obs.n_cols - 1);
              for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                tmpnm1(t) = alpha(i, t, k) + transition(i, j) + beta(j, t + 1, k);
                for (arma::uword r = 0; r < obs.n_rows; ++r) {
                  tmpnm1(t) += emission(j, obs(r, t + 1, k), r);
                }
              }
#pragma omp atomic
              ksii(i, j) += exp(logSumExp(tmpnm1) - ll(k));
            }
          }
        }
      }

      for (arma::uword r = 0; r < emission.n_slices; ++r) {
        for (arma::uword l = 0; l < nSymbols(r); ++l) {
          for (arma::uword i = 0; i < emission.n_rows; ++i) {
            if (emission(i, l, r) > -arma::datum::inf) {
              arma::vec tmpn(obs.n_cols);
              for (arma::uword t = 0; t < obs.n_cols; ++t) {
                if (l == (obs(r, t, k))) {
                  tmpn(t) = alpha(i, t, k) + beta(i, t, k);
                } else
                  tmpn(t) = -arma::datum::inf;
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
    for (arma::uword r = 0; r < emission.n_slices; ++r) {
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
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
    }

    double tmp = sum(ll);
    change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
    sumlogLik = tmp;
    if (!arma::is_finite(sumlogLik)) {
      return Rcpp::List::create(Rcpp::Named("error") = 6);
    }
    if (trace > 1) {
      Rcpp::Rcout << "iter: " << iter;
      Rcpp::Rcout << " logLik: " << sumlogLik;
      Rcpp::Rcout << " relative change: " << change << std::endl;
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
  return Rcpp::List::create(Rcpp::Named("initialProbs") = Rcpp::wrap(exp(init)),
      Rcpp::Named("transitionMatrix") = Rcpp::wrap(exp(transition)),
      Rcpp::Named("emissionArray") = Rcpp::wrap(exp(emission)), Rcpp::Named("logLik") = sumlogLik,
      Rcpp::Named("iterations") = iter, Rcpp::Named("change") = change, Rcpp::Named("error") = 0);
}

// EM algorithm for mixture hidden Markov models using logarithm space
// [[Rcpp::export]]
Rcpp::List log_EMx(const arma::mat& transition_, const arma::cube& emission_, 
                   const arma::vec& init_, const arma::ucube& obs, const arma::uvec& nSymbols, 
                   const arma::mat& coef_, const arma::mat& X, const arma::uvec& numberOfStates, 
                   int itermax, double tol, int trace, arma::uword threads) {
  
  // Make sure we don't alter the original vec/mat/cube
  // needed for cube, in future maybe in other cases as well
  arma::cube emission = log(emission_);
  arma::mat transition = log(transition_);
  arma::vec init = log(init_);
  arma::mat coef(coef_);
  
  coef.col(0).zeros();
  
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 3);
  }
  weights.each_row() /= sum(weights, 0);
  weights = log(weights);
  
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    initk.col(k) = init + reparma(weights.col(k), numberOfStates);
  }
  
  arma::cube alpha(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  arma::cube beta(emission.n_rows, obs.n_cols, obs.n_slices); //m,n,k
  
  log_internalForward(transition, emission, initk, obs, alpha, threads);
  log_internalBackward(transition, emission, obs, beta, threads);
  
  arma::vec ll(obs.n_slices);
  
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) num_threads(threads) \
  default(none) shared(obs, alpha, ll)
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
    }
    
    double sumlogLik = sum(ll);
  if (trace > 0) {
    Rcpp::Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
  }
  //  
  //  //EM-algorithm begins
  //  
  double change = tol + 1.0;
  int iter = 0;
  arma::uvec cumsumstate = cumsum(numberOfStates);
  
  while ((change > tol) && (iter < itermax)) {
    iter++;
    
    arma::mat ksii(emission.n_rows, emission.n_rows, arma::fill::zeros);
    arma::cube gamma(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
    arma::vec delta(emission.n_rows, arma::fill::zeros);
    
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      delta += exp(alpha.slice(k).col(0) + beta.slice(k).col(0) - ll(k));
    }
    
#pragma omp parallel for if(obs.n_slices>=threads) schedule(static) num_threads(threads) \
    default(none) shared(transition, obs, ll, alpha, beta, emission, ksii, gamma, nSymbols)
      for (arma::uword k = 0; k < obs.n_slices; ++k) {
        if (obs.n_cols > 1) {
          for (arma::uword j = 0; j < emission.n_rows; ++j) {
            for (arma::uword i = 0; i < emission.n_rows; ++i) {
              if (transition(i, j) > -arma::datum::inf) {
                arma::vec tmpnm1(obs.n_cols - 1);
                for (arma::uword t = 0; t < (obs.n_cols - 1); ++t) {
                  tmpnm1(t) = alpha(i, t, k) + transition(i, j) + beta(j, t + 1, k);
                  for (arma::uword r = 0; r < obs.n_rows; ++r) {
                    tmpnm1(t) += emission(j, obs(r, t + 1, k), r);
                  }
                }
#pragma omp atomic
                ksii(i, j) += exp(logSumExp(tmpnm1) - ll(k));
              }
            }
          }
        }
        
        for (arma::uword r = 0; r < emission.n_slices; ++r) {
          for (arma::uword l = 0; l < nSymbols(r); ++l) {
            for (arma::uword i = 0; i < emission.n_rows; ++i) {
              if (emission(i, l, r) > -arma::datum::inf) {
                arma::vec tmpn(obs.n_cols);
                for (arma::uword t = 0; t < obs.n_cols; ++t) {
                  if (l == (obs(r, t, k))) {
                    tmpn(t) = alpha(i, t, k) + beta(i, t, k);
                  } else
                    tmpn(t) = -arma::datum::inf;
                }
#pragma omp atomic
                gamma(i, l, r) += exp(logSumExp(tmpn) - ll(k));
                
              }
            }
          }
        }
      }
      
      arma::uword error = log_optCoef(weights, obs, emission, initk, beta, ll, coef, X, cumsumstate,
                                      numberOfStates, trace);
    if (error != 0) {
      return Rcpp::List::create(Rcpp::Named("error") = error);
    }
    if (obs.n_cols > 1) {
      ksii.each_col() /= sum(ksii, 1);
      transition = log(ksii);
    }
    for (arma::uword r = 0; r < emission.n_slices; ++r) {
      gamma.slice(r).cols(0, nSymbols(r) - 1).each_col() /= sum(
        gamma.slice(r).cols(0, nSymbols(r) - 1), 1);
      emission.slice(r).cols(0, nSymbols(r) - 1) = log(gamma.slice(r).cols(0, nSymbols(r) - 1));
    }
    
    for (arma::uword i = 0; i < numberOfStates.n_elem; ++i) {
      delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1) /= arma::as_scalar(
        arma::accu(delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1)));
    }
    
    init = log(delta);
    
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      initk.col(k) = init + reparma(weights.col(k), numberOfStates);
    }
    
    log_internalForward(transition, emission, initk, obs, alpha, threads);
    log_internalBackward(transition, emission, obs, beta, threads);
    
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      ll(k) = logSumExp(alpha.slice(k).col(obs.n_cols - 1));
    }
    
    double tmp = sum(ll);
    change = (tmp - sumlogLik) / (std::abs(sumlogLik) + 0.1);
    sumlogLik = tmp;
    if (!arma::is_finite(sumlogLik)) {
      return Rcpp::List::create(Rcpp::Named("error") = 6);
    }
    if (trace > 1) {
      Rcpp::Rcout << "iter: " << iter;
      Rcpp::Rcout << " logLik: " << sumlogLik;
      Rcpp::Rcout << " relative change: " << change << std::endl;
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
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::wrap(coef), Rcpp::Named("initialProbs") = Rcpp::wrap(exp(init)),
                            Rcpp::Named("transitionMatrix") = Rcpp::wrap(exp(transition)),
                            Rcpp::Named("emissionArray") = Rcpp::wrap(exp(emission)), Rcpp::Named("logLik") = sumlogLik,
                            Rcpp::Named("iterations") = iter, Rcpp::Named("change") = change, Rcpp::Named("error") = 0);
}
