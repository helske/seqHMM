#include "gradients.h"

arma::mat gradient_wrt_omega(
    const arma::mat& Qt, const arma::vec& omega, 
    const arma::vec& loglik_i, const arma::vec& loglik, const arma::mat& X, 
    const unsigned int i) {
  
  arma::mat gradmat = -omega * omega.t();
  gradmat.diag() += omega;
  arma::vec gradvec = gradmat * exp(loglik_i - loglik(i));
  return Qt * gradvec * X.col(i).t();
}

arma::mat gradient_wrt_pi(
    const arma::mat& Qt,
    const arma::mat& log_py, const arma::mat& log_beta, const double ll, 
    const arma::vec& Pi, const arma::mat& X, const unsigned int i) {
  
  unsigned int S = log_py.n_rows;
  arma::mat gradmat = -Pi * Pi.t();
  gradmat.diag() += Pi;
  arma::vec gradvec = gradmat * exp(log_py.col(0) + log_beta.col(0) - ll);
  
  return Qt * gradvec * X.col(i).t();
}
arma::mat gradient_wrt_pi(
    const arma::mat& Qt,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::vec>& Pi, const arma::mat& X, const unsigned int i,
    const unsigned int d) {
  
  unsigned int S = log_py.n_rows;
  arma::mat gradmat = -Pi(d) * Pi(d).t();
  gradmat.diag() += Pi(d);
  arma::vec gradvec = gradmat * exp(log_omega(d) + log_py.slice(d).col(0) + 
    log_beta.slice(d).col(0) - loglik(i));
  return Qt * gradvec * X.col(i).t();
}

arma::mat gradient_wrt_A(
    const arma::mat& Qt,
    const arma::mat& log_py, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& A, 
    const arma::cube& X, const unsigned int i, const unsigned int t, 
    const unsigned int s) {
  
  unsigned int S = log_py.n_rows;
  arma::mat gradmat = -A.slice(t).row(s).t() * A.slice(t).row(s);
  gradmat.diag() += A.slice(t).row(s);
  arma::vec gradvec = gradmat * exp(log_alpha(s, t) + log_py.col(t + 1) + log_beta.col(t + 1) - ll);
  return Qt * gradvec * X.slice(i).col(t).t();
}
arma::mat gradient_wrt_A(
    const arma::mat& Qt,
    const arma::vec& log_omega, const arma::cube& log_py, 
    const arma::cube& log_alpha, const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube> A,
    const arma::cube& X, const unsigned int i, const unsigned int t, 
    const unsigned int s, const unsigned int d) {
  
  unsigned int S = log_py.n_rows;
  arma::mat gradmat = -A(d).slice(t).row(s).t() * A(d).slice(t).row(s);
  gradmat.diag() += A(d).slice(t).row(s);
  arma::vec gradvec = gradmat * exp(log_omega(d) + log_alpha(s, t, d) + 
    log_py.slice(d).col(t + 1) + log_beta.slice(d).col(t + 1) - loglik(i));
  return Qt * gradvec * X.slice(i).col(t).t();
}
// NHMM singlechannel
arma::mat gradient_wrt_B_t0(
    const arma::mat& Qt,
    const arma::umat& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::cube& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s) {
  
  unsigned int M = B.n_cols - 1;
  arma::rowvec Brow = B.slice(0).row(s).cols(0, M - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(0, i));
  gradvec(obs(0, i)) += Brow(obs(0, i));
  double grad = exp(log_Pi(s) + log_beta(s, 0) - ll);
  return Qt * grad * gradvec * X.slice(i).col(0).t();
}
arma::mat gradient_wrt_B(
    const arma::mat& Qt,
    const arma::umat& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, 
    const arma::cube& log_A, const arma::cube& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s, const unsigned int t) {
  unsigned int M = B.n_cols - 1;
  arma::rowvec Brow = B.slice(t + 1).row(s).cols(0, M - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(t + 1, i));
  gradvec(obs(t + 1, i)) += Brow(obs(t + 1, i));
  double grad = arma::accu(
    exp(log_alpha.col(t) + log_A.slice(t).col(s) + log_beta(s, t + 1) - ll));
  return Qt * grad * gradvec * X.slice(i).col(t + 1).t();
}
// NHMM multichannel
arma::mat gradient_wrt_B_t0(
    const arma::mat& Qt,
    const arma::ucube& obs, const arma::vec& log_Pi, const arma::mat& log_beta, 
    const double ll, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const unsigned int i, const unsigned int s, 
    const unsigned int c) {
  
  unsigned int C = M.n_elem;
  arma::rowvec Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(c, 0, i));
  gradvec(obs(c, 0, i)) += Brow(obs(c, 0, i));
  double logpy = 0;
  for (unsigned int cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, 0, i), 0);
    }
  }
  double grad = exp(log_Pi(s) + logpy + log_beta(s, 0) - ll);
  return Qt * grad * gradvec * X.slice(i).col(0).t();
}
arma::mat gradient_wrt_B(
    const arma::mat& Qt,
    const arma::ucube& obs, const arma::mat& log_alpha, 
    const arma::mat& log_beta, const double ll, const arma::cube& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const unsigned int i, 
    const unsigned int s, const unsigned int t, const unsigned int c) {
  
  unsigned int C = M.n_elem;
  arma::rowvec Brow = B(c).slice(t + 1).row(s).cols(0, M(c) - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(c, t + 1, i));
  gradvec(obs(c, t + 1, i)) += Brow(obs(c, t + 1, i));
  double logpy = 0;
  for (unsigned int cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(cc, t + 1, i), t + 1);
    }
  }
  double grad = arma::accu(
    exp(log_alpha.col(t) + log_A.slice(t).col(s) + logpy + log_beta(s, t + 1) - ll));
  return Qt * grad * gradvec * X.slice(i).col(t + 1).t();
}
// MNHMM singlechannel
arma::mat gradient_wrt_B_t0(
    const arma::mat& Qt,
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& B, const arma::cube& X, 
    const unsigned int i, const unsigned int s, const unsigned d) {
  
  unsigned int M = B(d).n_cols - 1;
  arma::rowvec Brow = B(d).slice(0).row(s).cols(0, M - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(0, i));
  gradvec(obs(0, i)) += Brow(obs(0, i));
  double grad = exp(log_omega(d) + log_Pi(d)(s) + log_beta(s, 0, d) - loglik(i));
  return Qt * grad * gradvec * X.slice(i).col(0).t();
}
arma::mat gradient_wrt_B(
    const arma::mat& Qt,
    const arma::vec& log_omega,
    const arma::umat& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, const arma::field<arma::cube>& B, 
    const arma::cube& X, const unsigned int i, const unsigned int s, 
    const unsigned int t, const unsigned int d) {
  
  unsigned int M = B(0).n_cols - 1;
  arma::rowvec Brow = B(d).slice(t + 1).row(s).cols(0, M - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(t + 1, i));
  gradvec(obs(t + 1, i)) += Brow(obs(t + 1, i));
  double grad = arma::accu(
    exp(log_omega(d) + log_alpha.slice(d).col(t) + log_A(d).slice(t).col(s) + log_beta(s, t + 1, d) - loglik(i)));
  return Qt * grad * gradvec * X.slice(i).col(t + 1).t();
}
// MNHMM MC
arma::mat gradient_wrt_B_t0(
    const arma::mat& Qt,
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::field<arma::vec>& log_Pi, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, const arma::field<arma::cube>& log_B, 
    const arma::field<arma::cube>& B, const arma::cube& X, 
    const arma::uvec& M, const unsigned int i, const unsigned int s, 
    const unsigned int c, const unsigned int d) {
  
  unsigned int C = M.n_elem;
  arma::rowvec Brow = B(d * C + c).slice(0).row(s).cols(0, M(c) - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(c, 0, i));
  gradvec(obs(c, 0, i)) += Brow(obs(c, 0, i));
  double logpy = 0;
  for (unsigned int cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(d * C + cc)(s, obs(cc, 0, i), 0);
    }
  }
  double grad = exp(log_omega(d) + log_Pi(d)(s) + logpy + log_beta(s, 0, d) - loglik(i));
  return Qt * grad * gradvec * X.slice(i).col(0).t();
}

// MNHMM MC
arma::mat gradient_wrt_B(
    const arma::mat& Qt,
    const arma::vec& log_omega,
    const arma::ucube& obs, const arma::cube& log_alpha, 
    const arma::cube& log_beta, const arma::vec& loglik, 
    const arma::field<arma::cube>& log_A, 
    const arma::field<arma::cube>& log_B, const arma::field<arma::cube>& B, 
    const arma::cube& X, const arma::uvec& M, const unsigned int i, 
    const unsigned int s, const unsigned int t, const unsigned int c,
    const unsigned int d) {
  
  unsigned int C = M.n_elem;
  arma::rowvec Brow = B(d * C + c).slice(t + 1).row(s).cols(0, M(c) - 1);
  arma::vec gradvec = -Brow.t() * Brow(obs(c, t + 1, i));
  gradvec(obs(c, t + 1, i)) += Brow(obs(c, t + 1, i));
  double logpy = 0;
  for (unsigned int cc = 0; cc < C; cc++) {
    if (cc != c) {
      logpy += log_B(d * C + cc)(s, obs(cc, t + 1, i), t + 1);
    }
  }
  double grad = arma::accu(
    exp(log_omega(d) + log_alpha.slice(d).col(t) + 
      log_A(d).slice(t).col(s) + logpy + log_beta(s, t + 1, d) - loglik(i)));
  return Qt * grad * gradvec * X.slice(i).col(t + 1).t();
}