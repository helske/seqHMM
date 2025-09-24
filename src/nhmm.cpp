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
  py(S, Ti.max()),
  pi(S),
  A(S, S, Ti.max()),
  B(C),
  log_pi(S),
  log_A(S, S, Ti.max()),
  log_B(C),
  maxval(maxval),
  minval(
    (minval < 0)  ? std::pow(arma::datum::eps, 2.0/3.0) : minval
  ){
  for (arma::uword c = 0; c < C; ++c) {
    B(c) = arma::cube(S, M(c) + 1, Ti.max());
    log_B(c) = arma::cube(S, M(c) + 1, Ti.max());
  }
}

void nhmm::update_pi(const arma::uword i) {
  if (icpt_only_pi) {
    log_pi = log_softmax(gamma_pi.col(0));
    pi = softmax(gamma_pi.col(0));
  } else {
    arma::vec tmp = gamma_pi * X_pi.col(i);
    log_pi = log_softmax(tmp);
    pi = softmax(tmp);
  }
}
void nhmm::update_A(const arma::uword i) {
  arma::mat log_Atmp(S, S);
  arma::mat Atmp(S, S);
  if (icpt_only_A) {
    for (arma::uword s = 0; s < S; ++s) { // from states
      Atmp.col(s) = gamma_A.slice(s).col(0);
      log_Atmp.col(s) = log_softmax(Atmp.col(s));
      Atmp.col(s) = softmax(Atmp.col(s));
    }
    log_A.each_slice() = log_Atmp.t();
    A.each_slice() = Atmp.t();
  } else {
    if (tv_A) {
      for (arma::uword t = 0; t < Ti(i); ++t) { // time
        for (arma::uword s = 0; s < S; ++s) { // from states
          Atmp.col(s) = gamma_A.slice(s) * X_A(i).col(t);
          log_Atmp.col(s) = log_softmax(Atmp.col(s));
          Atmp.col(s) = softmax(Atmp.col(s));
        }
        log_A.slice(t) = log_Atmp.t();
        A.slice(t) = Atmp.t();
      }
    } else {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Atmp.col(s) = gamma_A.slice(s) * X_A(i).col(0);
        log_Atmp.col(s) = log_softmax(Atmp.col(s));
        Atmp.col(s) = softmax(Atmp.col(s));
      }
      log_A.each_slice() = log_Atmp.t();
      A.each_slice() = Atmp.t();
    }
  }
}
void nhmm::update_B(const arma::uword i) {
  for (arma::uword c = 0; c < C; ++c) {
    arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
    arma::mat log_Btmp(M(c) + 1, S, arma::fill::zeros);
    if (icpt_only_B(c)) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Btmp.col(s).rows(0, M(c) - 1) = softmax(
          gamma_B(c).slice(s).col(0)
        );
        log_Btmp.col(s).rows(0, M(c) - 1) = log_softmax(
          gamma_B(c).slice(s).col(0)
        );
      }
      log_B(c).each_slice() = log_Btmp.t();
      B(c).each_slice() = Btmp.t();
      
    } else {
      if (tv_B(c)) {
        for (arma::uword t = 0; t < Ti(i); ++t) { // time
          for (arma::uword s = 0; s < S; ++s) { // from states
            Btmp.col(s).rows(0, M(c) - 1) = gamma_B(c).slice(s) * X_B(c, i).col(t);
            log_Btmp.col(s).rows(0, M(c) - 1) = log_softmax(Btmp.col(s).rows(0, M(c) - 1));
            Btmp.col(s).rows(0, M(c) - 1) = softmax(Btmp.col(s).rows(0, M(c) - 1));
          }
          log_B(c).slice(t) = log_Btmp.t();
          B(c).slice(t) = Btmp.t();
        }
      } else {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s).rows(0, M(c) - 1) = gamma_B(c).slice(s) * X_B(c, i).col(0);
          log_Btmp.col(s).rows(0, M(c) - 1) = log_softmax(Btmp.col(s).rows(0, M(c) - 1));
          Btmp.col(s).rows(0, M(c) - 1) = softmax(Btmp.col(s).rows(0, M(c) - 1));
        }
        log_B(c).each_slice() = log_Btmp.t();
        B(c).each_slice() = Btmp.t();
      }
    }
  }
}
void nhmm::update_py(const arma::uword i) {
  py.ones();
  for (arma::uword t = 0; t < Ti(i); ++t) {
    for (arma::uword c = 0; c < C; ++c) {
      py.col(t) %= B(c).slice(t).col(obs(i)(c, t));
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
    update_py(i);
    logp(i) = univariate_viterbi(
      q(i), log_pi, log_A, arma::log(py.cols(0, Ti(i) - 1))
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
arma::vec nhmm::loglik() {
  arma::vec alpha(S);
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
    update_py(i);
    alpha = pi % py.col(0);
    double c = 1.0 / arma::accu(alpha);
    alpha *= c;
    ll(i) = -std::log(c);
    for (arma::uword t = 1; t < Ti(i); ++t) {
      alpha = A.slice(t).t() * alpha % py.col(t);
      c = 1.0 / arma::accu(alpha);
      alpha *= c;
      ll(i) -= std::log(c);
    }
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
    update_py(i);
    univariate_forward_logspace(
      log_alpha(i), log_pi, log_A, 
      arma::log(py.cols(0, Ti(i) - 1))
    );
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
    update_py(i);
    univariate_backward_logspace(
      log_beta(i), log_A, arma::log(py.cols(0, Ti(i) - 1))
    );
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
    arma::vec& tmpvec, 
    const arma::mat& beta,
    const arma::uword i) {
  tmpvec = py.col(0) % beta.col(0);
  grad += (pi % tmpvec - pi * arma::dot(pi, tmpvec)) * X_pi.col(i).t();
}

void nhmm::gradient_A(
    arma::mat& grad, 
    arma::vec& tmpvec1,
    arma::vec& tmpvec2,
    const arma::mat& alpha, 
    const arma::mat& beta, 
    const arma::uword i,
    const arma::uword t, 
    const arma::uword s) {
  
  tmpvec1 = A.slice(t).row(s).t();
  tmpvec2 = py.col(t) % beta.col(t) * alpha(s, t - 1);
  grad += (tmpvec1 %tmpvec2 - tmpvec1 * arma::dot(tmpvec1, tmpvec2)) * X_A(i).col(t).t();
}

void nhmm::gradient_B_t1(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::mat& beta, 
    const arma::uword i, 
    const arma::uword s, 
    const arma::uword c) {
  
  arma::subview_row<double> Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
  const arma::uword idx = obs(i)(c, 0);
  const double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double pyx = 1.0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      pyx *= B(cc)(s, obs(i)(cc, 0), 0);
    }
  }
  grad += pi(s) * pyx * beta(s, 0) * tmpvec * X_B(c, i).col(0).t();
}

void nhmm::gradient_B(
    arma::mat& grad, 
    arma::vec& tmpvec, 
    const arma::mat& alpha, 
    const arma::mat& beta, 
    const arma::uword i, 
    const arma::uword s,
    const arma::uword t,
    const arma::uword c) {
  
  arma::subview_row<double> Brow = B(c).slice(t).row(s).cols(0, M(c) - 1);
  const arma::uword idx = obs(i)(c, t);
  const double brow = Brow(idx);
  tmpvec = -Brow.t() * brow;
  tmpvec(idx) += brow;
  double pyx = 1.0;
  for (arma::uword cc = 0; cc < C; ++cc) {
    if (cc != c) {
      pyx *= B(cc)(s, obs(i)(cc, t), t);
    }
  }
  
  grad += arma::dot(alpha.col(t - 1), A.slice(t).col(s)) * beta(s, t) * 
    pyx * tmpvec * X_B(c, i).col(t).t();
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
  arma::vec tmpvecS1(S);
  arma::vec tmpvecS2(S);
  arma::field<arma::vec> tmpvecM(C);
  for (arma::uword c = 0; c < C; ++c) {
    tmpvecM(c) = arma::vec(M(c));
  }
  arma::vec loglik(N);
  arma::mat alpha(S, Ti.max());
  arma::mat beta(S, Ti.max());
  arma::vec scales(Ti.max());
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
    update_py(i);
    loglik(i) = univariate_forward(alpha, scales, pi, A, py, Ti(i));
    univariate_backward(beta, scales, A, py, Ti(i));
    if (!std::isfinite(loglik(i))) {
      return Rcpp::List::create(
        Rcpp::Named("loglik") = -maxval,
        Rcpp::Named("gradient_pi") = Rcpp::wrap(grad_pi2),
        Rcpp::Named("gradient_A") = Rcpp::wrap(grad_A2),
        Rcpp::Named("gradient_B") = Rcpp::wrap(grad_B2)
      );
    }
    // gradient wrt gamma_pi
    gradient_pi(grad_pi, tmpvecS1, beta, i);
    // gradient wrt gamma_A
    for (arma::uword t = 1; t < Ti(i); ++t) {
      for (arma::uword s = 0; s < S; ++s) {
        gradient_A(
          grad_A.slice(s), tmpvecS1, tmpvecS2, alpha, beta, i, t, s
        );
      }
    }
    for (arma::uword c = 0; c < C; ++c) {
      for (arma::uword s = 0; s < S; ++s) {
        if (obs(i)(c, 0) < M(c)) {
          gradient_B_t1(
            grad_B(c).slice(s), tmpvecM(c), beta, i, s, c
          );
        }
        for (arma::uword t = 1; t < Ti(i); ++t) {
          if (obs(i)(c, t) < M(c)) {
            gradient_B(
              grad_B(c).slice(s), tmpvecM(c), alpha, beta, i, s, t, c
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
