#include "config.h"
#include "mnhmm.h"
#include "softmax.h"
#include "eta_to_gamma.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"
#include "sample.h"

mnhmm::mnhmm(
  const arma::field<arma::umat>& obs,
  const arma::uvec& Ti,
  const arma::uvec& M,
  const arma::mat& X_pi,
  const arma::field<arma::mat>& X_A,
  arma::field<arma::mat>&& X_B,
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
  N(obs.n_elem),
  C(obs(0).n_rows),
  S(gamma_A(0).n_slices),
  D(gamma_A.n_elem),
  X_pi(X_pi),
  X_A(X_A),
  X_B(std::move(X_B)),
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
  log_py(S, Ti(0), D),
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
    A(d) = arma::cube(S, S, Ti(0));
    log_A(d) = arma::cube(S, S, Ti(0));
    for (arma::uword c = 0; c < C; ++c) {
      B(d, c) = arma::cube(S, M(c) + 1, Ti(0));
      log_B(d, c) = arma::cube(S, M(c) + 1, Ti(0));
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
  for (arma::uword d = 0; d < D; ++d) {
    A(d) = arma::cube(S, S, Ti(i));
  }
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
            Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A(i).col(t));
          }
          A(d).slice(t) = Atmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A(d).slice(s) * X_A(i).col(0));
        }
        A(d).each_slice() = Atmp.t();
      }
      log_A(d) = arma::log(A(d));
    }
  }
}
void mnhmm::update_B(const arma::uword i) {
  for (arma::uword c = 0; c < C; ++c) {
    arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
    for (arma::uword d = 0; d < D; ++d) {
      B(d, c) = arma::cube(S, M(c) + 1, Ti(i));
    }
    if (icpt_only_B(c)) {
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
      if (tv_B(c)) {
        for (arma::uword d = 0; d < D; ++d) {
          for (arma::uword t = 0; t < Ti(i); ++t) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(d, c).slice(s) * X_B(c, i).col(t));
            }
            B(d, c).slice(t) = Btmp.t();
          }
          log_B(d, c) = arma::log(B(d, c));
        }
      } else {
        for (arma::uword d = 0; d < D; ++d) {
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(
              gamma_B(d, c).slice(s) * X_B(c, i).col(0)
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
  log_py = arma::cube(S, Ti(i), D, arma::fill::zeros);
  for (arma::uword d = 0; d < D; ++d) {
    for (arma::uword t = 0; t < Ti(i); ++t) {
      for (arma::uword c = 0; c < C; ++c) {
        log_py.slice(d).col(t) += log_B(d, c).slice(t).col(obs(i)(c, t));
      }
    }
  }
}

Rcpp::List mnhmm::viterbi() {
  
  arma::field<arma::uvec> q(N);
  arma::vec logp(N, arma::fill::value(-arma::datum::inf));
  double logp_d;
  for (arma::uword i = 0; i < N; ++i) {
    q(i) = arma::uvec(Ti(i));
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
    arma::uvec q_d(Ti(i));
    for (arma::uword d = 0; d < D; ++d) {
      logp_d = univariate_viterbi(
        q_d, log_omega(d) + log_pi(d), log_A(d), log_py.slice(d)
      );
      if (logp_d > logp(i)) {
        logp(i) = logp_d;
        q(i) = q_d + d * S;
      }
    }
    
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}

arma::vec mnhmm::loglik() {
  arma::vec log_alpha(S);
  arma::vec log_alpha_new(S);
  arma::vec lli(D);
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
      log_alpha = log_pi(d) + log_py.slice(d).col(0);
      for (arma::uword t = 1; t < Ti(i); ++t) {
        for (arma::uword s = 0; s < S; ++s) {
          log_alpha_new(s) = logSumExp(
            log_alpha + log_A(d).slice(t).col(s) + log_py(s, t, d)
          );
        }
        log_alpha = log_alpha_new;
      }
      lli(d) = logSumExp(log_alpha);
    }
    ll(i) = logSumExp(log_omega + lli);
  }
  return ll;
}

arma::field<arma::mat> mnhmm::forward() {
  arma::field<arma::mat> log_alpha(N);
  for (arma::uword i = 0; i < N; ++i) {
    log_alpha(i) = arma::mat(D * S, Ti(i));
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
        log_alpha(i).rows(d * S, (d + 1) * S - 1);
      univariate_forward(
        submat, log_omega(d) + log_pi(d), log_A(d), log_py.slice(d)
      );
    }
  }
  return log_alpha;
}

arma::field<arma::mat> mnhmm::backward() {
  arma::field<arma::mat> log_beta(N);
  for (arma::uword i = 0; i < N; ++i) {
    log_beta(i) = arma::mat(D * S, Ti(i));
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    update_log_py(i);
    for (arma::uword d = 0; d < D; ++d) {
      arma::subview<double> submat = 
        log_beta(i).rows(d * S, (d + 1) * S - 1);
      univariate_backward(submat, log_A(d), log_py.slice(d));
    }
  }
  return log_beta;
}
Rcpp::List mnhmm::predict() {
  
  // these are P(z_t, y_t | data up to time t excluding y_t)
  arma::field<arma::cube> obs_prob(D, C, N);
  // these are P(z_t, z_{t-1} | data up to time t excluding y_t)
  arma::field<arma::cube> state_prob(D, N);
  for (arma::uword d = 0; d < D; ++d) {
    for (arma::uword i = 0; i < N; ++i) {
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(d, c, i) = arma::cube(S, M(c), Ti(i));
      }
      state_prob(d, i) = arma::cube(S, S, Ti(i));
    }
  }
  arma::vec alpha(S);
  arma::vec alpha_new(S);
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
      alpha = pi(d);
      alpha_new.ones();
      state_prob(d, i).slice(0).each_row() = alpha.t() / S;
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(d, c, i).slice(0) = 
          B(d, c).slice(0).cols(0, M(c) - 1).each_col() % alpha;
        if (obs(i)(c, 0) < M(c)) {
          alpha_new %= obs_prob(d, c, i).slice(0).col(obs(i)(c, 0));
        }
      }
      alpha = alpha_new / arma::accu(alpha_new);
      for (arma::uword t = 1; t < Ti(i); ++t) {
        state_prob(d, i).slice(t) = A(d).slice(t).each_col() % alpha;
        alpha = A(d).slice(t).t() * alpha ;
        alpha_new.ones();
        for (arma::uword c = 0; c < C; ++c) {
          obs_prob(d, c, i).slice(t) = 
            B(c).slice(t).cols(0, M(c) - 1).each_col() % alpha;
          if (obs(i)(c, t) < M(c)) {
            alpha_new %= obs_prob(d, c, i).slice(t).col(obs(i)(c, t));
          }
        }
        alpha = alpha_new / arma::accu(alpha_new);
      }
      // weight by cluster probability
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(d, c, i) *= omega(d); 
      }
      state_prob(d, i) *= omega(d);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(obs_prob),
    Rcpp::Named("states") = Rcpp::wrap(state_prob)
  );
}

Rcpp::List mnhmm::simulate() {
  
  arma::field<arma::umat> y(N);
  arma::field<arma::uvec> z(N);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::uvec seqD = arma::linspace<arma::uvec>(0, D - 1, D);
  arma::field<arma::uvec> seqM(C);
  for (arma::uword c = 0; c < C; ++c) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  for (arma::uword i = 0; i < N; ++i) {
    y(i) = arma::umat(C, Ti(i));
    z(i) = arma::uvec(Ti(i));
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
    z(i)(0) = sample(seqS, pi(cluster).t());
    for (arma::uword c = 0; c < C; ++c) {
      y(i)(c, 0) = sample(seqM(c), B(cluster, c).slice(0).row(z(i)(0)).cols(0, M(c) - 1));
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      z(i)(t) = sample(seqS, A(cluster).slice(t).row(z(i)(t - 1)));
      for (arma::uword c = 0; c < C; ++c) {
        y(i)(c, t) = sample(seqM(c), B(cluster, c).slice(t).row(z(i)(t)).cols(0, M(c) - 1));
      }
    }
    z(i) += cluster * S;
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
    loglik(i)) * X_A(i).col(t).t();
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
  arma::uword idx = obs(i)(c, 0);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(d, cc)(s, obs(i)(cc, 0), 0);
    }
  }
  grad += exp(log_omega(d) + log_pi(d)(s) + logpy + 
    log_beta(s, 0, d) - loglik(i)) * tmpvec * X_B(c, i).col(0).t();
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
  arma::uword idx = obs(i)(c, t);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(d, cc)(s, obs(i)(cc, t), t);
    }
  }
  grad += arma::accu(exp(log_omega(d) + log_alpha.slice(d).col(t - 1) + 
    log_A(d).slice(t).col(s) + logpy + log_beta(s, t, d) - loglik(i))) * 
    tmpvec * X_B(c, i).col(t).t();
}

Rcpp::List mnhmm::log_objective(const arma::mat& Qs, 
                                const arma::field<arma::mat>& Qm,
                                const arma::mat& Qd) {
  
  arma::mat grad_omega(D, X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi(D);
  arma::field<arma::cube> grad_A(D);
  arma::field<arma::cube> grad_B(C, D);
  for (arma::uword d = 0; d < D; ++d) {
    grad_pi(d) = arma::mat(S, X_pi.n_rows, arma::fill::zeros);
    grad_A(d) = arma::cube(S, X_A(0).n_rows, S, arma::fill::zeros);
    for (arma::uword c = 0; c < C; ++c) {
      grad_B(c, d) = arma::cube(M(c), X_B(c, 0).n_rows, S, arma::fill::zeros);
    }
  }
  arma::mat grad_omega2(D - 1, X_omega.n_rows, arma::fill::zeros);
  arma::field<arma::mat> grad_pi2(D);
  arma::field<arma::cube> grad_A2(D);
  arma::field<arma::cube> grad_B2(C, D);
  arma::vec loglik(N);
  arma::vec loglik_i(D);
  arma::cube log_alpha(S, Ti(0), D);
  arma::cube log_beta(S, Ti(0), D);
  for (arma::uword d = 0; d < D; ++d) {
    grad_pi2(d) = arma::mat(S - 1, X_pi.n_rows, arma::fill::zeros);
    grad_A2(d) = arma::cube(S - 1, X_A(0).n_rows, S, arma::fill::zeros);
    for (arma::uword c = 0; c < C; ++c) {
      grad_B2(c, d) = arma::cube(M(c) - 1, X_B(c, 0).n_rows, S, arma::fill::zeros);
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
    log_alpha = arma::cube(S, Ti(i), D);
    log_beta = arma::cube(S, Ti(i), D);
    for (arma::uword d = 0; d < D; ++d) {
      univariate_forward(
        log_alpha.slice(d), log_pi(d), log_A(d), log_py.slice(d)
      );
      univariate_backward(log_beta.slice(d), log_A(d), log_py.slice(d));
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
          if (obs(i)(c, 0) < M(c)) {
            gradient_B_t1(
              grad_B(c, d).slice(s), tmpvec(c), log_beta, loglik, i, s, c, d
            );
          }
          for (arma::uword t = 1; t < Ti(i); ++t) {
            if (obs(i)(c, t) < M(c)) {
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

