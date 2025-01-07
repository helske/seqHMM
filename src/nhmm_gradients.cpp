#include "nhmm_gradients.h"

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
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& pi, const arma::mat& X, const arma::uword i) {
  
  tmpmat = -pi * pi.t();
  tmpmat.diag() += pi;
  grad += tmpmat * exp(log_py.col(0) + log_beta.col(0) - ll) * X.col(i).t();
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
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const arma::uword i, const arma::uword t, 
    const arma::uword s) {
  
  tmpmat = -A.slice(t).row(s).t() * A.slice(t).row(s);
  tmpmat.diag() += A.slice(t).row(s);
  grad += tmpmat * exp(log_alpha(s, t) + log_py.col(t + 1) + 
    log_beta.col(t + 1) - ll) * X.slice(i).col(t).t();
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
  grad += tmpmat * exp(log_omega(d) + log_alpha(s, t, d) + 
    log_py.slice(d).col(t + 1) + log_beta.slice(d).col(t + 1) - 
    loglik(i)) * X.slice(i).col(t).t();
}
// NHMM singlechannel
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::umat& obs, const arma::vec& log_pi, const arma::mat& log_beta, 
    const double ll, const arma::cube& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s) {
  
  arma::rowvec Brow = B.slice(0).row(s).cols(0, B.n_cols - 2);
  arma::uword idx = obs(0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  grad += exp(log_pi(s) + log_beta(s, 0) - ll) * tmpvec * X.slice(i).col(0).t();
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::umat& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, 
    const arma::cube& log_A, const arma::cube& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s, const arma::uword t) {
  
  arma::rowvec Brow = B.slice(t + 1).row(s).cols(0, B.n_cols - 2);
  arma::uword idx = obs(t + 1, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  grad += arma::accu(exp(log_alpha.col(t) + log_A.slice(t).col(s) + 
    log_beta(s, t + 1) - ll)) * tmpvec * X.slice(i).col(t + 1).t();
}
// NHMM multichannel
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::vec& log_pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, 0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, 0, i), 0);
    }
  }
  grad += exp(log_pi(s) + logpy + log_beta(s, 0) - ll) * tmpvec * 
    X.slice(i).col(0).t();
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(t + 1).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, t + 1, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, t + 1, i), t + 1);
    }
  }
  grad += arma::accu(exp(log_alpha.col(t) + log_A.slice(t).col(s) + 
    logpy + log_beta(s, t + 1) - ll)) * tmpvec * X.slice(i).col(t + 1).t();
}
// MNHMM singlechannel
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::field<arma::vec>& log_pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uword i, const arma::uword s, const unsigned d) {
  
  arma::uword M = B(d).n_cols - 1;
  arma::rowvec Brow = B(d).slice(0).row(s).cols(0, M - 1);
  arma::uword idx = obs(0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  grad += exp(log_omega(d) + log_pi(d)(s) + log_beta(s, 0, d) - 
    loglik(i)) * tmpvec * X.slice(i).col(0).t();
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uword i, const arma::uword s, 
    const arma::uword t, const arma::uword d) {
  
  arma::uword M = B(0).n_cols - 1;
  arma::rowvec Brow = B(d).slice(t + 1).row(s).cols(0, M - 1);
  arma::uword idx = obs(t + 1, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  grad += arma::accu(exp(log_omega(d) + log_alpha.slice(d).col(t) + 
    log_A(d).slice(t).col(s) + log_beta(s, t + 1, d) - loglik(i))) * 
    tmpvec * X.slice(i).col(t + 1).t();
}
// MNHMM MC
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const arma::uword i, const arma::uword s, 
    const arma::uword c, const arma::uword d) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(d * C + c).slice(0).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, 0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(d * C + cc)(s, obs(cc, 0, i), 0);
    }
  }
  grad += exp(log_omega(d) + log_pi(d)(s) + logpy + 
    log_beta(s, 0, d) - loglik(i)) * tmpvec * X.slice(i).col(0).t();
}

// MNHMM MC
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, 
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const arma::uword i, 
    const arma::uword s, const arma::uword t, const arma::uword c,
    const arma::uword d) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(d * C + c).slice(t + 1).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(c, t + 1, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(d * C + cc)(s, obs(cc, t + 1, i), t + 1);
    }
  }
  grad += arma::accu(exp(log_omega(d) + log_alpha.slice(d).col(t) + 
    log_A(d).slice(t).col(s) + logpy + log_beta(s, t + 1, d) - loglik(i))) * 
    tmpvec * X.slice(i).col(t + 1).t();
}

//FANHMM
void gradient_wrt_A(
    arma::mat& grad, arma::mat& tmpmat,  const arma::umat& obs, 
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const arma::cube& W, const arma::uword i, 
    const arma::uword t, const arma::uword s) {
  arma::uword n = W.n_rows; 
  tmpmat = -A.slice(t).row(s).t() * A.slice(t).row(s);
  tmpmat.diag() += A.slice(t).row(s);
  arma::vec tmpvec = tmpmat * exp(log_alpha(s, t) + log_py.col(t + 1) + 
    log_beta.col(t + 1) - ll);
  grad.cols(0, X.n_rows - 1) += tmpvec * X.slice(i).col(t).t();
  if (obs(t, i) > 0 && n > 0) {
    grad.cols(X.n_rows + (obs(t, i) - 1) * n, X.n_rows + obs(t, i) * n - 1) += 
      tmpvec * W.slice(i).col(t).t();
  }
}

// FANHMM singlechannel
void gradient_wrt_B_t0(
    arma::mat& grad, arma::vec& tmpvec, const arma::umat& obs, 
    const arma::uvec& obs_0, const arma::vec& log_pi, const arma::mat& log_beta, 
    const double ll, const arma::cube& B, const arma::cube& X, 
    const arma::cube& W, const arma::uword i, const arma::uword s) {
  
  arma::uword n = W.n_rows; 
  arma::rowvec Brow = B.slice(0).row(s).cols(0, B.n_cols - 2);
  arma::uword idx = obs(0, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  tmpvec = exp(log_pi(s) + log_beta(s, 0) - ll) * tmpvec;
  grad.cols(0, X.n_rows - 1) += tmpvec * X.slice(i).col(0).t();
  if (obs_0(i) > 0 && n > 0) {
    grad.cols(X.n_rows + (obs_0(i) - 1) * n, X.n_rows + obs_0(i) * n - 1) += 
      tmpvec * W.slice(i).col(0).t();
  }
}
void gradient_wrt_B(
    arma::mat& grad, arma::vec& tmpvec, const arma::umat& obs, 
    const arma::mat& log_alpha, const arma::mat& log_beta, const double ll, 
    const arma::cube& log_A, const arma::cube& B, const arma::cube& X, 
    const arma::cube& W, const arma::uword i, const arma::uword s, 
    const arma::uword t) {
  
  arma::uword n = W.n_rows; 
  arma::rowvec Brow = B.slice(t + 1).row(s).cols(0, B.n_cols - 2);
  arma::uword idx = obs(t + 1, i);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  tmpvec = arma::accu(exp(log_alpha.col(t) + log_A.slice(t).col(s) + 
    log_beta(s, t + 1) - ll)) * tmpvec;
  grad.cols(0, X.n_rows - 1) += tmpvec * X.slice(i).col(t + 1).t();
  if (obs(t, i) > 0 && n > 0) {
    grad.cols(X.n_rows + (obs(t, i) - 1) * n, X.n_rows + obs(t, i) * n - 1) += 
      tmpvec * W.slice(i).col(t + 1).t();
  }
}