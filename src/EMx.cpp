// EM algorithm for mixture hidden Markov models
#include "forward_backward.h"
#include "optcoef.h"
#include "reparma.h"
#include <omp.h>

// [[Rcpp::export]]
Rcpp::List EMx(const arma::mat& transition_, const arma::cube& emission_, const arma::vec& init_,
  const arma::ucube& obs, const arma::uvec& nSymbols, const arma::mat& coef_, const arma::mat& X,
  const arma::uvec& numberOfStates, int itermax, double tol, int trace, unsigned int threads) {
  
  // Make sure we don't alter the original vec/mat/cube
  // needed for cube, in future maybe in other cases as well
  arma::cube emission(emission_);
  arma::mat transition(transition_);
  arma::vec init(init_);
  arma::mat coef(coef_);
  
  coef.col(0).zeros();
  arma::mat weights = exp(X * coef).t();
  if (!weights.is_finite()) {
    return Rcpp::List::create(Rcpp::Named("error") = 3);
  }
  weights.each_row() /= sum(weights, 0);
  
  arma::mat initk(emission.n_rows, obs.n_slices);
  for (unsigned int k = 0; k < obs.n_slices; k++) {
    initk.col(k) = init % reparma(weights.col(k), numberOfStates);
  }
  
  //
  //  //EM-algorithm begins
  //
  double change = tol + 1.0;
  int iter = 0;
  
  arma::uvec cumsumstate = arma::cumsum(numberOfStates);
  double sumlogLik_new = 0;
  double sumlogLik = -1e150;
  
  while ((change > tol) && (iter < itermax)) {
    iter++;
    
    arma::mat ksii(emission.n_rows, emission.n_rows, arma::fill::zeros);
    arma::cube gamma(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
    arma::vec delta(emission.n_rows, arma::fill::zeros);
    
    arma::mat bsi(emission.n_rows, obs.n_slices);
    sumlogLik_new = 0;
    double max_sf = 1;
    unsigned int error_code = 0;
    
#pragma omp parallel for if(obs.n_slices >= threads) schedule(static) reduction(+:sumlogLik_new) num_threads(threads) \
    default(shared) //shared(bsi, initk, transition, obs, emission, delta, ksii, gamma, nSymbols, error_code, max_sf, arma::fill::zeros)
      for (unsigned int k = 0; k < obs.n_slices; k++) {
        
        if (error_code == 0) {
          arma::mat alpha(emission.n_rows, obs.n_cols); //m,n,k
          arma::vec scales(obs.n_cols);
          arma::sp_mat sp_trans(transition);
          uvForward(sp_trans.t(), emission, initk.col(k), obs.slice(k), alpha, scales);
          arma::mat beta(emission.n_rows, obs.n_cols); //m,n,k
          uvBackward(sp_trans, emission, obs.slice(k), beta, scales);
          sumlogLik_new -= arma::sum(log(scales));
          
          arma::mat ksii_k(emission.n_rows, emission.n_rows, arma::fill::zeros);
          arma::cube gamma_k(emission.n_rows, emission.n_cols, emission.n_slices, arma::fill::zeros);
          arma::vec delta_k(emission.n_rows);
          delta_k = alpha.col(0) % beta.col(0) / scales(0);
          
          for (unsigned int i = 0; i < emission.n_rows; i++) {
            for (unsigned int j = 0; j < emission.n_rows; j++) {
              if (transition(i, j) > 0.0) {
                for (unsigned int t = 0; t < (obs.n_cols - 1); t++) {
                  double tmp = alpha(i, t) * transition(i, j) * beta(j, t + 1);
                  for (unsigned int r = 0; r < obs.n_rows; r++) {
                    tmp *= emission(j, obs(r, t + 1, k), r);
                  }
                  ksii_k(i, j) += tmp;
                }
              }
            }
          }
          for (unsigned int r = 0; r < emission.n_slices; r++) {
            for (unsigned int l = 0; l < nSymbols(r); l++) {
              for (unsigned int i = 0; i < emission.n_rows; i++) {
                if (emission(i, l, r) > 0.0) {
                  for (unsigned int t = 0; t < obs.n_cols; t++) {
                    if (l == (obs(r, t, k))) {
                      double tmp = alpha(i, t) * beta(i, t) / scales(t);
                      gamma_k(i, l, r) += tmp;
                    }
                  }
                }
              }
            }
          }
          
          for (unsigned int j = 0; j < emission.n_rows; j++) {
            bsi(j, k) = beta(j, 0)  * initk(j, k);
          }
          
#pragma omp critical
{
  if(!scales.is_finite()) {
    error_code = 1;
  }
  if(!beta.is_finite()) {
    error_code = 2;
  }
  max_sf = std::min(max_sf, scales.max());
  delta += delta_k;
  ksii += ksii_k;
  gamma += gamma_k;
}
        }
      }
      if(error_code == 1) {
        return Rcpp::List::create(Rcpp::Named("error") = 1);
      }
      if(error_code == 2) {
        return Rcpp::List::create(Rcpp::Named("error") = 2);
      }
      if (max_sf > 1e150) {
        Rcpp::warning("Largest scaling factor was %e, results can be numerically unstable.", max_sf);
      }
      change = (sumlogLik_new - sumlogLik) / (std::abs(sumlogLik) + 0.1);
      sumlogLik = sumlogLik_new;
      
      if (trace > 0) {
        if(iter == 0) {
          Rcpp::Rcout << "Log-likelihood of initial model: " << sumlogLik << std::endl;
        } else {
          if (trace > 1) {
            Rcpp::Rcout << "iter: " << iter;
            Rcpp::Rcout << " logLik: " << sumlogLik;
            Rcpp::Rcout << " relative change: " << change << std::endl;
          }
        }
      }
      if (change > tol) {
        unsigned int error = optCoef(weights, obs, emission, bsi, coef, X, cumsumstate,
          numberOfStates, trace);
        if (error != 0) {
          return Rcpp::List::create(Rcpp::Named("error") = error);
        }
        
        if (obs.n_cols > 1) {
          arma::vec rsums = sum(ksii, 1);
          for (unsigned int kk = 0; kk < ksii.n_rows; kk++) {
            if (rsums(kk) == 0) {
              rsums(kk) = 1;
            }
          }
          ksii.each_col() /= rsums;
          transition = ksii;
        }
        for (unsigned int r = 0; r < emission.n_slices; r++) {
          
          gamma.slice(r).cols(0, nSymbols(r) - 1).each_col() /= sum(
            gamma.slice(r).cols(0, nSymbols(r) - 1), 1);
          emission.slice(r).cols(0, nSymbols(r) - 1) = gamma.slice(r).cols(0, nSymbols(r) - 1);
        }
        
        for (unsigned int i = 0; i < numberOfStates.n_elem; i++) {
          delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1) /= arma::as_scalar(
            arma::accu(delta.subvec(cumsumstate(i) - numberOfStates(i), cumsumstate(i) - 1)));
        }
        init = delta;
        
        for (unsigned int k = 0; k < obs.n_slices; k++) {
          initk.col(k) = init % reparma(weights.col(k), numberOfStates);
        }
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
  return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::wrap(coef), Rcpp::Named("initialProbs") = Rcpp::wrap(init),
    Rcpp::Named("transitionMatrix") = Rcpp::wrap(transition), Rcpp::Named("emissionArray") = Rcpp::wrap(emission),
    Rcpp::Named("logLik") = sumlogLik, Rcpp::Named("iterations") = iter, Rcpp::Named("change") = change,
      Rcpp::Named("error") = 0);
}
