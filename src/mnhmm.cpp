#include "config.h"
#include "mnhmm.h"
#include "softmax.h"
#include "eta_to_gamma.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"
#include "sample.h"

mnhmm::mnhmm(
  const arma::ucube& obs,
  const arma::uvec& Ti,
  const arma::uvec& M,
  const arma::mat& X_pi,
  const arma::cube& X_A,
  const arma::field<arma::cube>& X_B,
  const arma::mat& X_omega,
  const bool icpt_only_pi,
  const bool icpt_only_A,
  const arma::uvec& icpt_only_B,
  const bool icpt_only_omega,
  const bool iv_A,
  const arma::uvec& iv_B,
  const bool tv_A,
  const arma::uvec& tv_B,
  const arma::field<arma::mat>& gamma_pi,
  const arma::field<arma::cube>& gamma_A,
  const arma::field<arma::cube>& gamma_B,
  const arma::mat& gamma_omega,
  double maxval,
  double minval)
  :
  obs(obs),
  Ti(Ti),
  M(M),
  N(obs.n_slices),
  T(obs.n_cols),
  C(obs.n_rows),
  S(gamma_A(0).n_slices),
  D(gamma_A.n_elem),
  X_pi(X_pi),
  X_A(X_A),
  X_B(X_B),
  X_omega(X_omega),
  icpt_only_pi(icpt_only_pi),
  icpt_only_A(icpt_only_A),
  icpt_only_B(icpt_only_B),
  icpt_only_omega(icpt_only_omega),
  iv_A(iv_A),
  iv_B(iv_B),
  tv_A(tv_A),
  tv_B(tv_B),
  gamma_pi(gamma_pi),
  gamma_A(gamma_A),
  gamma_B(gamma_B),
  gamma_omega(gamma_omega),
  log_py(S, T, D),
  pi(D),
  log_pi(D),
  A(D),
  log_A(D),
  B(D, C),
  log_B(D, C),
  omega(D),
  log_omega(D),
  maxval(maxval),
  minval(
    (minval < 0)  ? std::pow(arma::datum::eps, 2.0/3.0) : minval
  )
{
  for (arma::uword d = 0; d < D; ++d) {
    pi(d) = arma::vec(S);
    log_pi(d) = arma::vec(S);
    A(d) = arma::cube(S, S, T);
    log_A(d) = arma::cube(S, S, T);
    for (arma::uword c = 0; c < C; ++c) {
      B(d, c) = arma::cube(S, M(c) + 1, T);
      log_B(d, c) = arma::cube(S, M(c) + 1, T);
    }
  }
}

void mnhmm::update_pi(const arma::uword i) {
  if (icpt_only_pi) {
    for (arma::uword d = 0; d < D; ++d) {
      pi(d) = softmax(gamma_pi(d).col(0));
      log_pi(d) = arma::log(pi(d));
    }
  } else {
    for (arma::uword d = 0; d < D; ++d) {
      pi(d) = softmax(gamma_pi(d) * X_pi.col(i));
      log_pi(d) = arma::log(pi(d));
    }
  }
}
void mnhmm::update_A(const arma::uword i) {
  arma::mat Atmp(S, S);
  if (icpt_only_A) {
    for (arma::uword d = 0; d < D; ++d) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Atmp.col(s) = softmax(gamma_A(d).slice(s).col(0));
      }
      A(d).each_slice() = Atmp.t();
      log_A(d) = arma::log(A(d));
    }
  } else {
    for (arma::uword d = 0; d < D; ++d) {
      if (tv_A) {
        for (arma::uword t = 0; t < Ti(i); ++t) { // time
          for (arma::uword s = 0; s < S; ++s) { // from states
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(t));
          }
          A(d).slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A.slice(i).col(0));
        }
        A(d).each_slice() = Atmp.t();
      }
      log_A(d) = arma::log(A(d));
    }
  }
}
void mnhmm::update_B(const arma::uword i) {
  for (arma::uword c = 0; c < C; ++c) {
    if (icpt_only_B(c)) {
      arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
      for (arma::uword d = 0; d < D; ++d) {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(d, c).slice(s).col(0)
          );
        }
        B(d, c).each_slice() = Btmp.t();
        log_B(d, c) = arma::log(B(d, c));
      }
    } else { 
      arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
      if (tv_B(c)) {
        for (arma::uword d = 0; d < D; ++d) {
          for (arma::uword t = 0; t < Ti(i); ++t) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(d, c).slice(s) * X_B(c).slice(i).col(t));
            }
            B(d, c).slice(t) = Btmp.t();
          }
          log_B(d, c) = arma::log(B(d, c));
        }
      } else {
        for (arma::uword d = 0; d < D; ++d) {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(d, c).slice(s) * X_B(c).slice(i).col(0)
            );
          }
          B(d, c).each_slice() = Btmp.t();
          log_B(d, c) = arma::log(B(d, c));
        }
      }
    }
  }
}
void mnhmm::update_omega(const arma::uword i) {
  if (icpt_only_omega) {
    omega = softmax(gamma_omega.col(0));
  } else {
    omega = softmax(gamma_omega * X_omega.col(i));  
  }
  log_omega = arma::log(omega);
}
void mnhmm::update_log_py(const arma::uword i) {
  log_py.zeros();
  for (arma::uword d = 0; d < D; ++d) {
    for (arma::uword t = 0; t < Ti(i); ++t) {
      for (arma::uword c = 0; c < C; ++c) {
        log_py.slice(d).col(t) += log_B(d, c).slice(t).col(obs(c, t, i));
      }
    }
  }
}

void mnhmm::viterbi(arma::umat& q, arma::vec& logp) {
  logp.fill(-arma::datum::inf);
  double logp_d;
  arma::uvec q_d(q.n_rows);
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      logp_d = univariate_viterbi(
        q_d,
        log_omega(d) + log_pi(d),
        log_A(d).slices(0, Ti(i) - 1),
        log_py.slice(d).cols(0, Ti(i) - 1)
      );
      if (logp_d > logp(i)) {
        logp(i) = logp_d;
        q.col(i) = q_d;
      }
    }
  }
}
arma::vec mnhmm::loglik() {
  arma::vec log_alpha(S);
  arma::vec log_alpha_new(S);
  arma::mat lli(S, D);
  arma::vec ll(N, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      log_alpha = log_omega(d) + log_pi(d) + log_py.slice(d).col(0);
      for (arma::uword t = 1; t < Ti(i); ++t) {
        for (arma::uword s = 0; i < S; ++s) {
          log_alpha_new(s) = logSumExp(
            log_alpha + log_A(d).slice(t).col(s) + log_py(s, t, d)
          );
        }
        log_alpha = log_alpha_new;
      }
      lli.col(d) = logSumExp(log_alpha);
    }
    ll(i) = logSumExp(arma::vectorise(lli));
  }
  return ll;
}

arma::cube mnhmm::forward() {
  arma::cube log_alpha(
      S * D, T, N, arma::fill::value(arma::datum::nan)
  );
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      arma::subview<double> submat = 
        log_alpha.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_forward(
        submat,
        log_omega(d) + log_pi(d),
        log_A(d), 
        log_py.slice(d).cols(0, Ti(i) - 1)
      );
    }
  }
  return log_alpha;
}

arma::cube mnhmm::backward() {
  arma::cube log_beta(
      S * D, T, N, arma::fill::value(arma::datum::nan)
  );
  for (arma::uword i = 0; i < N; ++i) {
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      arma::subview<double> submat = 
        log_beta.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_backward(
        submat,
        log_A(d), 
        log_py.slice(d).cols(0, Ti(i) - 1)
      );
    }
  }
  return log_beta;
}
Rcpp::List mnhmm::predict() {
  
  // these are P(z_t, y_t | data up to time t excluding y_t)
  arma::field<arma::cube> obs_prob(C, N);
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword i = 0; i < N; ++i) {
      obs_prob(c, i) = arma::cube(S * D, M(c), Ti(i));
    }
  }
  // these are P(z_t, z_{t-1} | data up to time t excluding y_t)
  arma::field<arma::cube> state_prob(N);
  for (arma::uword i = 0; i < N; ++i) {
    state_prob(i) = arma::cube(S * D, S * D, Ti(i));
  }
  arma::vec alpha(S);
  arma::vec alpha_new(S);
  arma::field<arma::span> blockM(C);
  for (arma::uword c = 0; c < C; ++c) {
    blockM(c) = arma::span(0, M(c) - 1);
  }
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    for (arma::uword d = 0; d < D; ++d) {
      arma::span block(d * S, (d + 1) * S - 1);
      alpha = pi(d);
      alpha_new.ones();
      state_prob(i).slice(0)(block, block).each_row() = alpha.t() / S;
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(c, i).slice(0)(block, blockM(c)) = 
          B(c).slice(0).cols(blockM(c)).each_col() % alpha;
        if (obs(c, 0, i) < M(c)) {
          alpha_new %= obs_prob(c, i).slice(0).col(obs(c, 0, i));
        }
      }
      alpha = alpha_new / arma::accu(alpha_new);
      for (arma::uword t = 1; t < Ti(i); ++t) {
        state_prob(i).slice(t)(block, block) = A(d).slice(t).each_col() % alpha;
        alpha = A(d).slice(t).t() * alpha ;
        alpha_new.ones();
        for (arma::uword c = 0; c < C; ++c) {
          obs_prob(c, i).slice(t)(block, blockM(c)) = 
            B(c).slice(t).cols(blockM(c)).each_col() % alpha;
          if (obs(c, t, i) < M(c)) {
            alpha_new %= obs_prob(c, i).slice(t).col(obs(c, t, i));
          }
        }
        alpha = alpha_new / arma::accu(alpha_new);
      }
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(c, i).rows(d * S, (d + 1) * S - 1) *= omega(d); 
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(obs_prob),
    Rcpp::Named("states") = Rcpp::wrap(state_prob)
  );
}

Rcpp::List mnhmm::simulate() {
  
  arma::ucube y(C, T, N, arma::fill::value(arma::datum::nan));
  arma::umat z(T, N, arma::fill::value(arma::datum::nan));
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  arma::field<arma::uvec> seqM(C);
  for (arma::uword c = 0; c < C; ++c) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  for (arma::uword i = 0; i < N; ++i) {
    
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    arma::uword cluster = sample(seqD, omega.t());
    z(0, i) = sample(seqS, pi(cluster).t());
    for (arma::uword c = 0; c < C; ++c) {
      y(c, 0, i) = sample(seqM(c), B(cluster, c).slice(0).row(z(0, i)).cols(0, M(c) - 1));
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      z(t, i) = sample(seqS, A(cluster).slice(t).row(z(t - 1, i)));
      for (arma::uword c = 0; c < C; ++c) {
        y(c, t, i) = sample(seqM(c), B(cluster, c).slice(t).row(z(t, i)).cols(0, M(c) - 1));
      }
    }
    z.col(i) += cluster * S;
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y),
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
void mnhmm::gradient_omega(
    arma::mat& grad, arma::mat& tmpmat,
    const arma::vec& loglik_i, const arma::vec& loglik,
    const arma::uword i) {
  
  tmpmat = -omega * omega.t();
  tmpmat.diag() += omega;
  grad += tmpmat * exp(loglik_i - loglik(i)) * X_omega.col(i).t();
}
void mnhmm::gradient_pi(
    arma::mat& grad, 
    arma::mat& tmpmat, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, 
    const arma::uword i,
    const arma::uword d) {
  
  tmpmat = -pi(d) * pi(d).t();
  tmpmat.diag() += pi(d);
  grad += tmpmat * exp(log_omega(d) + log_py.slice(d).col(0) + 
    log_beta.slice(d).col(0) - loglik(i)) * X_pi.col(i).t();
}
void mnhmm::gradient_A(
    arma::mat& grad, 
    arma::mat& tmpmat,
    const arma::cube& log_alpha, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, 
    const arma::uword i, 
    const arma::uword t, 
    const arma::uword s, 
    const arma::uword d) {
  
  tmpmat = -A(d).slice(t).row(s).t() * A(d).slice(t).row(s);
  tmpmat.diag() += A(d).slice(t).row(s);
  grad += tmpmat * exp(log_omega(d) + log_alpha(s, t - 1, d) + 
    log_py.slice(d).col(t) + log_beta.slice(d).col(t) - 
    loglik(i)) * X_A.slice(i).col(t).t();
}
void mnhmm::gradient_B_t1(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, 
    const arma::uword i, 
    const arma::uword s, 
    const arma::uword c, 
    const arma::uword d) {
  
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
    log_beta(s, 0, d) - loglik(i)) * tmpvec * X_B(c).slice(i).col(0).t();
}
void mnhmm::gradient_B(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::cube& log_alpha, 
    const arma::cube& log_beta, 
    const arma::vec& loglik, 
    const arma::uword i, 
    const arma::uword s, 
    const arma::uword t, 
    const arma::uword c,
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
    tmpvec * X_B(c).slice(i).col(t).t();
}

Rcpp::List mnhmm::log_objective(const arma::mat& Qs, 
                                const arma::field<arma::mat>& Qm,
                                const arma::mat& Qd) {
  
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, T, D);
  arma::cube log_beta(S, T, D);
  arma::mat grad_omega(D, X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(D);
  arma::field<arma::cube> grad_A(D);
  arma::field<arma::cube> grad_B(C, D);
  for (arma::uword d = 0; d < D; ++d) {
    grad_pi(d) = arma::mat(S, X_pi.n_rows, arma::fill::zeros);
    grad_A(d) = arma::cube(S, X_A.n_rows, S, arma::fill::zeros);
    for (arma::uword c = 0; c < C; ++c) {
      grad_B(c, d) = arma::cube(M(c), X_B(c).n_rows, S, arma::fill::zeros);
    }
  }
  arma::mat grad_omega2(D - 1, X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi2(D);
  arma::field<arma::cube> grad_A2(D);
  arma::field<arma::cube> grad_B2(C, D);
  for (arma::uword d = 0; d < D; ++d) {
    grad_pi2(d) = arma::mat(S - 1, X_pi.n_rows, arma::fill::zeros);
    grad_A2(d) = arma::cube(S - 1, X_A.n_rows, S, arma::fill::zeros);
    for (arma::uword c = 0; c < C; ++c) {
      grad_B2(c, d) = arma::cube(M(c) - 1, X_B(c).n_rows, S, arma::fill::zeros);
    }
  }
  arma::mat tmpmat(S, S);
  arma::mat tmpmatD(D, D);
  arma::field<arma::vec> tmpvec(C);
  for (arma::uword c = 0; c < C; ++c) {
    tmpvec(c) = arma::vec(M(c));
  }
  for (arma::uword i = 0; i < N; ++i) {
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    if (!icpt_only_omega || i == 0) {
      update_omega(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      univariate_forward(
        log_alpha.slice(d), log_pi(d), log_A(d), 
        log_py.slice(d).cols(0, Ti(i) - 1)
      );
      univariate_backward(
        log_beta.slice(d), log_A(d), log_py.slice(d).cols(0, Ti(i) - 1)
      );
      loglik_i(d) = logSumExp(log_alpha.slice(d).col(Ti(i) - 1));
    }
    loglik(i) = logSumExp(log_omega + loglik_i);
    if (!std::isfinite(loglik(i))) {
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -maxval,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2),
        Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega2)
      );
    }
    // gradient wrt gamma_pi
    // d loglik / d pi
    for (arma::uword d = 0; d < D; ++d) {
      // gradient wrt gamma_pi
      gradient_pi(grad_pi(d), tmpmat, log_beta, loglik, i, d );
      // gradient wrt gamma_A
      for (arma::uword t = 1; t < Ti(i); ++t) {
        for (arma::uword s = 0; s < S; ++s) {
          gradient_A(
            grad_A(d).slice(s), tmpmat, log_alpha, log_beta, loglik, i, t, s, d
          );
        }
      }
      // gradient wrt gamma_B
      for (arma::uword c = 0; c < C; ++c) {
        for (arma::uword s = 0; s < S; ++s) {
          if (obs(c, 0, i) < M(c)) {
            gradient_B_t1(
              grad_B(c, d).slice(s), tmpvec(c), log_beta, loglik, i, s, c, d
            );
          }
          for (arma::uword t = 1; t < Ti(i); ++t) {
            if (obs(c, t, i) < M(c)) {
              gradient_B(
                grad_B(c, d).slice(s), tmpvec(c), log_alpha, log_beta, loglik, 
                i, s, t, c, d
              );
            }
          }
        }
      }
    }
    gradient_omega(
      grad_omega, tmpmatD, loglik_i, loglik, i
    );
  }
  grad_omega2 = Qd.t() * grad_omega;
  for (arma::uword d = 0; d < D; ++d) {
    grad_pi2(d) = Qs.t() * grad_pi(d);
    for (arma::uword s = 0; s < S; ++s) {
      grad_A2(d).slice(s) = Qs.t() * grad_A(d).slice(s);
    }
  }
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword d = 0; d < D; ++d) {
      for (arma::uword s = 0; s < S; ++s) {
        grad_B2(c, d).slice(s) =  Qm(c).t() * grad_B(c, d).slice(s);
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2),
    Rcpp::Named("gradient_omega") = Rcpp::wrap(grad_omega2)
  );
}

