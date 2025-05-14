#include "config.h"
#include "nhmm.h"
#include "softmax.h"
#include "eta_to_gamma.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"
#include "sample.h"

nhmm::nhmm(
  const arma::field<arma::umat>& obs,
  const arma::uvec& Ti,
  const arma::uvec& M,
  const arma::mat& X_pi,
  const arma::field<arma::mat>& X_A,
  arma::field<arma::mat>&& X_B,
  const bool icpt_only_pi,
  const bool icpt_only_A,
  const arma::uvec& icpt_only_B,
  const bool iv_A,
  const arma::uvec& iv_B,
  const bool tv_A,
  const arma::uvec& tv_B,
  const arma::mat& gamma_pi,
  const arma::cube& gamma_A,
  const arma::field<arma::cube>& gamma_B,
  double maxval, double minval)
  :
  obs(obs),
  Ti(Ti),
  M(M),
  N(obs.n_elem),
  C(obs(0).n_rows),
  S(gamma_A.n_slices),
  X_pi(X_pi),
  X_A(X_A),
  X_B(std::move(X_B)),
  icpt_only_pi(icpt_only_pi),
  icpt_only_A(icpt_only_A),
  icpt_only_B(icpt_only_B),
  iv_A(iv_A),
  iv_B(iv_B),
  tv_A(tv_A),
  tv_B(tv_B),
  gamma_pi(gamma_pi),
  gamma_A(gamma_A),
  gamma_B(gamma_B),
  log_py(S, Ti(0)),
  pi(S),
  log_pi(S),
  A(S, S, Ti(0)),
  log_A(S, S, Ti(0)),
  B(C),
  log_B(C),
  maxval(maxval),
  minval(
    (minval < 0)  ? std::pow(arma::datum::eps, 2.0/3.0) : minval
  ){
  for (arma::uword c = 0; c < C; ++c) {
    B(c) = arma::cube(S, M(c) + 1, Ti(0));
    log_B(c) = arma::cube(S, M(c) + 1, Ti(0));
  }
}

void nhmm::update_pi(const arma::uword i) {
  if (icpt_only_pi) {
    pi = softmax(gamma_pi.col(0));
  } else {
    pi = softmax(gamma_pi * X_pi.col(i));
  }
  log_pi = arma::log(pi);
}
void nhmm::update_A(const arma::uword i) {
  A = arma::cube(S, S, Ti(i));
  arma::mat Atmp(S, S);
  if (icpt_only_A) {
    for (arma::uword s = 0; s < S; ++s) { // from states
      Atmp.col(s) = softmax(gamma_A.slice(s).col(0));
    }
    A.each_slice() = Atmp.t();
  } else {
    if (tv_A) {
      for (arma::uword t = 0; t < Ti(i); ++t) { // time
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = softmax(gamma_A.slice(s) * X_A(i).col(t));
        }
        A.slice(t) = Atmp.t();
      }
    } else {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Atmp.col(s) = softmax(gamma_A.slice(s) * X_A(i).col(0));
      }
      A.each_slice() = Atmp.t();
    }
  }
  log_A = arma::log(A);
}
void nhmm::update_B(const arma::uword i) {
  for (arma::uword c = 0; c < C; ++c) {
    B(c) = arma::cube(S, M(c) + 1, Ti(i));
    arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
    if (icpt_only_B(c)) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Btmp.col(s).rows(0, M(c) - 1) = softmax(
          gamma_B(c).slice(s).col(0)
        );
      }
      B(c).each_slice() = Btmp.t();
    } else {
      if (tv_B(c)) {
        for (arma::uword t = 0; t < Ti(i); ++t) { // time
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = softmax(gamma_B(c).slice(s) * X_B(c, i).col(t));
          }
          B(c).slice(t) = Btmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = softmax(
            gamma_B(c).slice(s) * X_B(c, i).col(0)
          );
        }
        B(c).each_slice() = Btmp.t();
      }
    }
    log_B(c) = arma::log(B(c));
  }
}
void nhmm::update_log_py(const arma::uword i) {
  log_py = arma::mat(S, Ti(i), arma::fill::zeros);
  for (arma::uword t = 0; t < Ti(i); ++t) {
    for (arma::uword c = 0; c < C; ++c) {
      log_py.col(t) += log_B(c).slice(t).col(obs(i)(c, t));
    }
  }
}

Rcpp::List nhmm::viterbi() {
  arma::field<arma::uvec> q(N);
  arma::vec logp(N);
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
    update_log_py(i);
    logp(i) = univariate_viterbi(q(i), log_pi, log_A, log_py);
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
arma::vec nhmm::loglik() {
  arma::vec log_alpha(S);
  arma::vec log_alpha_new(S);
  arma::vec ll(N);
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
    update_log_py(i);
    log_alpha = log_pi + log_py.col(0);
    for (arma::uword t = 1; t < Ti(i); ++t) {
      for (arma::uword s = 0; s < S; ++s) {
        log_alpha_new(s) = logSumExp(
          log_alpha + log_A.slice(t).col(s) + log_py(s, t)
        );
      }
      log_alpha = log_alpha_new;
    }
    ll(i) = logSumExp(log_alpha);
  }
  return ll;
}

arma::field<arma::mat> nhmm::forward() {
  arma::field<arma::mat> log_alpha(N);
  for (arma::uword i = 0; i < N; ++i) {
    log_alpha(i) = arma::mat(S, Ti(i));
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    update_log_py(i);
    univariate_forward(log_alpha(i), log_pi, log_A, log_py);
  }
  return log_alpha;
}

arma::field<arma::mat>  nhmm::backward() {
  arma::field<arma::mat> log_beta(N);
  for (arma::uword i = 0; i < N; ++i) {
    log_beta(i) = arma::mat(S, Ti(i));
    if (!icpt_only_pi || i == 0) {
      update_pi(i);
    }
    if (iv_A || i == 0) {
      update_A(i);
    }
    if (arma::any(iv_B) || i == 0) {
      update_B(i);
    }
    update_log_py(i);
    univariate_backward(log_beta(i),log_A, log_py);
  }
  return log_beta;
}

Rcpp::List nhmm::predict() {
  
  // these are P(z_t, y_t | data up to time t excluding y_t)
  arma::field<arma::cube> obs_prob(C, N);
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword i = 0; i < N; ++i) {
      obs_prob(c, i) = arma::cube(S, M(c), Ti(i));
    }
  }
  // these are P(z_t, z_{t-1} | data up to time t excluding y_t)
  arma::field<arma::cube> state_prob(N);
  for (arma::uword i = 0; i < N; ++i) {
    state_prob(i) = arma::cube(S, S, Ti(i));
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
    alpha = pi;
    alpha_new.ones();
    state_prob(i).slice(0).each_row() = alpha.t() / S;
    for (arma::uword c = 0; c < C; ++c) {
      obs_prob(c, i).slice(0) = B(c).slice(0).cols(0, M(c) - 1).each_col() % alpha;
      if (obs(i)(c, 0) < M(c)) {
        alpha_new %= obs_prob(c, i).slice(0).col(obs(i)(c, 0));
      }
    }
    alpha = alpha_new / arma::accu(alpha_new);
    for (arma::uword t = 1; t < Ti(i); ++t) {
      state_prob(i).slice(t) = A.slice(t).each_col() % alpha;
      alpha = A.slice(t).t() * alpha;
      alpha_new.ones();
      for (arma::uword c = 0; c < C; ++c) {
        obs_prob(c, i).slice(t) = B(c).slice(t).cols(0, M(c) - 1).each_col() % alpha;
        if (obs(i)(c, t) < M(c)) {
          alpha_new %= obs_prob(c, i).slice(t).col(obs(i)(c, t));
        }
      }
      alpha = alpha_new / arma::accu(alpha_new);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(obs_prob),
    Rcpp::Named("states") = Rcpp::wrap(state_prob)
  );
}


Rcpp::List nhmm::simulate() {
  
  arma::field<arma::umat> y(N);
  arma::field<arma::uvec> z(N);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
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
    z(i)(0) = sample(seqS, pi.t());
    for (arma::uword c = 0; c < C; ++c) {
      y(i)(c, 0) = sample(seqM(c), B(c).slice(0).row(z(i)(0)).cols(0, M(c) - 1));
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      z(i)(t) = sample(seqS, A.slice(t).row(z(i)(t - 1)));
      for (arma::uword c = 0; c < C; ++c) {
        y(i)(c, t) = sample(seqM(c), B(c).slice(t).row(z(i)(t)).cols(0, M(c) - 1));
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y),
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
void nhmm::gradient_pi(
    arma::mat& grad, 
    arma::mat& tmpmat, 
    const arma::mat& log_beta,
    const double ll,
    const arma::uword i) {
  
  tmpmat = -pi * pi.t();
  tmpmat.diag() += pi;
  grad += tmpmat * exp(log_py.col(0) + log_beta.col(0) - ll) * X_pi.col(i).t();
}

void nhmm::gradient_A(
    arma::mat& grad, 
    arma::mat& tmpmat,
    const arma::mat& log_alpha, 
    const arma::mat& log_beta, 
    const double ll,  
    const arma::uword i,
    const arma::uword t, 
    const arma::uword s) {
  
  tmpmat = -A.slice(t).row(s).t() * A.slice(t).row(s);
  tmpmat.diag() += A.slice(t).row(s);
  grad += tmpmat * exp(log_alpha(s, t - 1) + log_py.col(t) + 
    log_beta.col(t) - ll) * X_A(i).col(t).t();
}

void nhmm::gradient_B_t1(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::mat& log_beta, 
    const double ll, 
    const arma::uword i, 
    const arma::uword s, 
    const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(i)(c, 0);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(i)(cc, 0), 0);
    }
  }
  grad += exp(log_pi(s) + logpy + log_beta(s, 0) - ll) * tmpvec * 
    X_B(c, i).col(0).t();
}
void nhmm::gradient_B(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::mat& log_alpha, 
    const arma::mat& log_beta, 
    const double ll,
    const arma::uword i, 
    const arma::uword s,
    const arma::uword t,
    const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow = B(c).slice(t).row(s).cols(0, M(c) - 1);
  arma::uword idx = obs(i)(c, t);
  double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double logpy = 0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      logpy += log_B(cc)(s, obs(i)(cc, t), t);
    }
  }
  grad += arma::accu(exp(log_alpha.col(t - 1) + log_A.slice(t).col(s) + 
    logpy + log_beta(s, t) - ll)) * tmpvec * X_B(c, i).col(t).t();
}
Rcpp::List nhmm::log_objective(const arma::mat& Qs, 
                               const arma::field<arma::mat>& Qm) {
 
  arma::mat grad_pi(S, X_pi.n_rows, arma::fill::zeros);
  arma::cube grad_A(S, X_A(0).n_rows, S, arma::fill::zeros);
  arma::field<arma::cube> grad_B(C);
  for (arma::uword c = 0; c < C; ++c) {
    grad_B(c) = arma::cube(M(c), X_B(c, 0).n_rows, S, arma::fill::zeros);
  }
  arma::mat grad_pi2(S - 1, X_pi.n_rows, arma::fill::zeros);
  arma::cube grad_A2(S - 1, X_A(0).n_rows, S, arma::fill::zeros);
  arma::field<arma::cube> grad_B2(C);
  for (arma::uword c = 0; c < C; ++c) {
    grad_B2(c) = arma::cube(M(c) - 1, X_B(c, 0).n_rows, S, arma::fill::zeros);
  }
  arma::mat tmpmat(S, S);
  arma::field<arma::vec> tmpvec(C);
  for (arma::uword c = 0; c < C; ++c) {
    tmpvec(c) = arma::vec(M(c));
  }
  arma::vec loglik(N);
  arma::mat log_alpha(S, Ti(0));
  arma::mat log_beta(S, Ti(0));
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
    update_log_py(i);
    log_alpha = arma::mat(S, Ti(i));
    log_beta = arma::mat(S, Ti(i));
    univariate_forward(log_alpha, log_pi, log_A, log_py);
    univariate_backward(log_beta, log_A, log_py);
    double ll = logSumExp(log_alpha.col(Ti(i) - 1));
    if (!std::isfinite(ll)) {
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -maxval,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2)
      );
    }
    loglik(i) = ll;
    // gradient wrt gamma_pi
    gradient_pi(grad_pi, tmpmat, log_beta, ll, i);
    // gradient wrt gamma_A
    for (arma::uword t = 1; t < Ti(i); ++t) {
      for (arma::uword s = 0; s < S; ++s) {
        gradient_A(
          grad_A.slice(s), tmpmat, log_alpha, log_beta, ll, i, t, s
        );
      }
    }
   
    for (arma::uword c = 0; c < C; ++c) {
      for (arma::uword s = 0; s < S; ++s) {
        if (obs(i)(c, 0) < M(c)) {
          gradient_B_t1(
            grad_B(c).slice(s), tmpvec(c), log_beta, ll, i, s, c
          );
        }
        for (arma::uword t = 1; t < Ti(i); ++t) {
          if (obs(i)(c, t) < M(c)) {
            gradient_B(
              grad_B(c).slice(s), tmpvec(c), log_alpha, log_beta,
              ll, i, s, t, c
            );
          }
        }
      }
    }
  }
  grad_pi2 = Qs.t() * grad_pi;
  for (arma::uword s = 0; s < S; ++s) {
    grad_A2.slice(s) = Qs.t() * grad_A.slice(s);
  }
  for (arma::uword c = 0; c < C; ++c) {
    for (arma::uword s = 0; s < S; ++s) {
      grad_B2(c).slice(s) = Qm(c).t() * grad_B(c).slice(s);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("loglik") = sum(loglik),
    Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
    Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
    Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2)
  );
}
