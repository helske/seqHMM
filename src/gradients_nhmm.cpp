#include "gradients_nhmm.h"

void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat, 
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& pi, const arma::mat& X, const arma::uword i) {
  
  tmpmat = -pi * pi.t();
  tmpmat.diag() += pi;
  grad += tmpmat * exp(log_py.col(0) + log_beta.col(0) - ll) * X.col(i).t();
}

void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s) {
  
  tmpmat = -A.slice(t).row(s).t() * A.slice(t).row(s);
  tmpmat.diag() += A.slice(t).row(s);
  grad += tmpmat * exp(log_alpha(s, t - 1) + log_py.col(t) + 
    log_beta.col(t) - ll) * X.slice(i).col(t).t();
}

void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::vec& log_pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::field<arma::cube>& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, 0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, 0, i), 0);
    }
  }
  grad += exp(log_pi(s) + logpy + log_beta(s, 0) - ll) * tmpvec * 
    X(c).slice(i).col(0).t();
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::field<arma::cube>& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(t).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, t, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, t, i), t);
    }
  }
  grad += arma::accu(exp(log_alpha.col(t - 1) + log_A.slice(t).col(s) + 
    logpy + log_beta(s, t) - ll)) * tmpvec * X(c).slice(i).col(t).t();
}
