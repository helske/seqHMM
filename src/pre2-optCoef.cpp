//Estimation of beta coefficients
#include "pre2-optcoef.h"

//Variance-covariance matrix of gamma coefficients based on inverse of -Hessian
// [[Rcpp::export]]
Rcpp::NumericMatrix varcoef(const arma::mat& coef, const arma::mat& X) {
  
  arma::mat weights = exp(X * coef).t();
  weights.each_row() /= sum(weights, 0);
  // use inv instead of faster inv_sympd as the latter produces error on valgrind
  return Rcpp::wrap(arma::inv(-hCoef(weights, X)));
}

arma::uword optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission,
    const arma::mat& bsi, arma::mat& coef,
    const arma::mat& X, const arma::uvec& cumsumstate, const arma::uvec& numberOfStates,
    int trace) {

  int iter = 0;
  double change = 1.0;
  while ((change > 1e-10) && (iter < 100)) {
    arma::vec tmpvec(X.n_cols * (weights.n_rows - 1));
    bool solve_ok = arma::solve(tmpvec, hCoef(weights, X),
        gCoef(obs, bsi, emission, weights, X, cumsumstate, numberOfStates),
        arma::solve_opts::no_approx);
    if (solve_ok == false) {
      return (4);
    }

    arma::mat coefnew(coef.n_rows, coef.n_cols - 1);
    for (arma::uword i = 0; i < (weights.n_rows - 1); ++i) {
      coefnew.col(i) = coef.col(i + 1) - tmpvec.subvec(i * X.n_cols, (i + 1) * X.n_cols - 1);
    }
    change = arma::accu(arma::abs(coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) - coefnew))
        / coefnew.n_elem;
    coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) = coefnew;
    iter++;
    if (trace == 3) {
      Rcpp::Rcout << "coefficient optimization iter: " << iter;
      Rcpp::Rcout << " new coefficients: " << std::endl << coefnew << std::endl;
      Rcpp::Rcout << " relative change: " << change << std::endl;
    }
    weights = exp(X * coef).t();
    if (!weights.is_finite()) {
      return (5);
    }
    weights.each_row() /= sum(weights, 0);

  }
  return (0);
}

arma::vec gCoef(const arma::ucube& obs, const arma::mat& bsi,
    const arma::cube& emission, const arma::mat& weights,
    const arma::mat& X, const arma::uvec& cumsumstate, const arma::uvec& numberOfStates) {

  int q = X.n_cols;
  arma::vec grad(q * (weights.n_rows - 1), arma::fill::zeros);
  double tmp;

  for (arma::uword k = 0; k < obs.n_slices; ++k) {
    for (arma::uword jj = 1; jj < numberOfStates.n_elem; ++jj) {
      for (arma::uword j = 0; j < emission.n_rows; ++j) {
        tmp = 1.0;
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          tmp *= emission(j, obs(r, 0, k), r);
        }
        if ((j >= (cumsumstate(jj) - numberOfStates(jj))) && (j < cumsumstate(jj))) {
          grad.subvec(q * (jj - 1), q * jj - 1) += tmp * bsi(j, k) *
               X.row(k).t() * (1.0 - weights(jj, k));
        } else {
          grad.subvec(q * (jj - 1), q * jj - 1) -= tmp * bsi(j,k) * X.row(k).t() * weights(jj, k);
        }
      }
    }
  }

  return grad;
}

arma::mat hCoef(const arma::mat& weights, const arma::mat& X) {

  int p = X.n_cols;
  arma::mat hess(p * (weights.n_rows - 1), p * (weights.n_rows - 1), arma::fill::zeros);
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    arma::mat XX = X.row(i).t() * X.row(i);
    for (arma::uword j = 0; j < (weights.n_rows - 1); ++j) {
      for (arma::uword k = 0; k < (weights.n_rows - 1); ++k) {
        if (j != k) {
          hess.submat(j * p, k * p, (j + 1) * p - 1, (k + 1) * p - 1) += XX * weights(j + 1, i)
              * weights(k + 1, i);
        } else {
          hess.submat(j * p, j * p, (j + 1) * p - 1, (j + 1) * p - 1) -= XX * weights(j + 1, i)
              * (1.0 - weights(j + 1, i));

        }

      }
    }
  }
  return hess;
}

// estimation of gamma coefficients using log-space

arma::uword log_optCoef(arma::mat& weights, const arma::ucube& obs, const arma::cube& emission,
                        const arma::mat& initk, const arma::cube& beta, const arma::vec& ll, arma::mat& coef,
                        const arma::mat& X, const arma::uvec& cumsumstate, const arma::uvec& numberOfStates,
                        int trace) {
  
  weights = exp(X * coef).t();
  weights.each_row() /= sum(weights, 0);
  int iter = 0;
  double change = 1.0;
  while ((change > 1e-10) && (iter < 100)) {
    arma::vec tmpvec(X.n_cols * (weights.n_rows - 1));
    bool solve_ok = arma::solve(tmpvec, hCoef(weights, X),
                                log_gCoef(obs, beta, emission, initk, weights, ll, X, cumsumstate, numberOfStates),
                                arma::solve_opts::no_approx);
    
    if (solve_ok == false) {
      return (4);
    }
    
    arma::mat coefnew(coef.n_rows, coef.n_cols - 1);
    for (arma::uword i = 0; i < (weights.n_rows - 1); ++i) {
      coefnew.col(i) = coef.col(i + 1) - tmpvec.subvec(i * X.n_cols, (i + 1) * X.n_cols - 1);
    }
    change = arma::accu(arma::abs(coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) - coefnew))
      / coefnew.n_elem;
    coef.submat(0, 1, coef.n_rows - 1, coef.n_cols - 1) = coefnew;
    iter++;
    if (trace == 3) {
      Rcpp::Rcout << "coefficient optimization iter: " << iter;
      Rcpp::Rcout << " new coefficients: " << std::endl << coefnew << std::endl;
      Rcpp::Rcout << " relative change: " << change << std::endl;
    }
    weights = exp(X * coef).t();
    if (!weights.is_finite()) {
      return (5);
    }
    weights.each_row() /= sum(weights, 0);
    
  }
  weights = log(weights);
  return (0);
}

arma::vec log_gCoef(const arma::ucube& obs, const arma::cube& beta, const arma::cube& emission,
                    const arma::mat& initk, const arma::mat& weights, const arma::vec& ll, const arma::mat& X,
                    const arma::uvec& cumsumstate, const arma::uvec& numberOfStates) {
  
  int q = X.n_cols;
  arma::vec grad(q * (weights.n_rows - 1), arma::fill::zeros);
  double tmp;
  for (arma::uword jj = 1; jj < numberOfStates.n_elem; ++jj) {
    for (arma::uword k = 0; k < obs.n_slices; ++k) {
      for (arma::uword j = 0; j < emission.n_rows; ++j) {
        tmp = 0.0;
        for (arma::uword r = 0; r < obs.n_rows; ++r) {
          tmp += emission(j, obs(r, 0, k), r);
        }
        if ((j >= (cumsumstate(jj) - numberOfStates(jj))) && (j < cumsumstate(jj))) {
          grad.subvec(q * (jj - 1), q * jj - 1) += exp(tmp + beta(j, 0, k) - ll(k) + initk(j, k))
          * X.row(k).t() * (1.0 - weights(jj, k));
        } else {
          grad.subvec(q * (jj - 1), q * jj - 1) -= exp(tmp + beta(j, 0, k) - ll(k) + initk(j, k))
          * X.row(k).t() * weights(jj, k);
        }
        
      }
    }
  }
  
  return grad;
}
