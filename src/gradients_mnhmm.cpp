#include "gradients_mnhmm.h"

void gradient_wrt_omega(
    arma::mat& grad, arma::mat& tmpmat, 
    const arma::vec& omega, 
    const arma::vec& loglik_i, const arma::vec& loglik, const arma::mat& X, 
    const arma::uword i) {
  
  tmpmat = -omega * omega.t();
  tmpmat.diag() += omega;
  grad += tmpmat * exp(loglik_i - loglik(i)) * X.col(i).t();
}
void gradient_wrt_pi(
    arma::mat& grad, arma::mat& tmpmat, 
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::vec>& pi, const arma::mat& X, const arma::uword i,
    const arma::uword d) {
  
  tmpmat = -pi(d) * pi(d).t();
  tmpmat.diag() += pi(d);
  grad += tmpmat * exp(log_omega(d) + log_py.slice(d).col(0) + 
    log_beta.slice(d).col(0) - loglik(i)) * X.col(i).t();
}
void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat, 
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_alpha, const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube> A,
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s, const arma::uword d) {
  
  tmpmat = -A(d).slice(t).row(s).t() * A(d).slice(t).row(s);
  tmpmat.diag() += A(d).slice(t).row(s);
  grad += tmpmat * exp(log_omega(d) + log_alpha(s, t - 1, d) + 
    log_py.slice(d).col(t) + log_beta.slice(d).col(t) - 
    loglik(i)) * X.slice(i).col(t).t();
}
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::field<arma::cube>& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c, const arma::uword d) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(d, c).slice(0).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, 0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(d, cc)(s, obs(cc, 0, i), 0);
    }
  }
  grad += exp(log_omega(d) + log_pi(d)(s) + logpy + 
    log_beta(s, 0, d) - loglik(i)) * tmpvec * X(c).slice(i).col(0).t();
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::field<arma::cube>& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c,
    const arma::uword d) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(d, c).slice(t).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, t, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(d, cc)(s, obs(cc, t, i), t);
    }
  }
  grad += arma::accu(exp(log_omega(d) + log_alpha.slice(d).col(t - 1) + 
    log_A(d).slice(t).col(s) + logpy + log_beta(s, t, d) - loglik(i))) * 
    tmpvec * X(c).slice(i).col(t).t();
}
