#include "config.h"
#include "fanhmm.h"
#include "softmax.h"
#include "list_to_field.h"
#include "joint_probability.h"
#include "sample.h"
#include "get_W_idx.h"

fanhmm::fanhmm(
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
  const arma::vec& prior_y,
  const Rcpp::List& W_X_B,
  double maxval,
  double minval)
  : nhmm(obs, Ti, M, X_pi, X_A, std::move(X_B), icpt_only_pi, icpt_only_A, 
    icpt_only_B, iv_A, iv_B, tv_A, tv_B, gamma_pi, gamma_A, gamma_B, 
    maxval, minval), prior_y(prior_y), fixed_0(prior_y.n_elem == 1),
    W_X_B(fixed_0 ? arma::field<arma::vec>() : veclist_to_3d_field(W_X_B)),
    B1(C), log_B1(C) {
  for (arma::uword c = 0; c < C; ++c) {
    B1(c) = arma::cube(S, M(c) + 1, prior_y.n_elem);
    log_B1(c) = arma::cube(S, M(c) + 1, prior_y.n_elem);
  }
}
void fanhmm::update_B(const arma::uword i) {
  for (arma::uword c = 0; c < C; ++c) {
    arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
    B(c) = arma::cube(S, M(c) + 1, Ti(i));
    if (icpt_only_B(c)) {
      for (arma::uword s = 0; s < S; ++s) { // from states
        Btmp.col(s).rows(0, M(c) - 1) = softmax(
          gamma_B(c).slice(s).col(0)
        );
      }
      B(c).each_slice() = Btmp.t();
      log_B(c) = arma::log(B(c));
    } else {
      if (tv_B(c)) {
        if (!fixed_0) {
          B1(c).zeros();
          B(c).slice(0).cols(0, M(c) - 1).zeros();
          for (arma::uword j = 0; j < prior_y.n_elem; ++j) {
            for (arma::uword s = 0; s < S; ++s) { // from states
              B1(c).slice(j).row(s).cols(0, M(c) - 1) = softmax(
                gamma_B(c).slice(s) * W_X_B(j, c, i)
              ).t();
            }
            B(c).slice(0) += B1(c).slice(j) * prior_y(j);
            B1(c).col(M(c)).ones();
          }
          log_B1(c) = arma::log(B1(c));
        }
        for (arma::uword t = 1 - fixed_0; t < Ti(i); ++t) { // time
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
      log_B(c) = arma::log(B(c));
    }
  }
}

Rcpp::List fanhmm::predict(
    arma::field<arma::mat>&& W_A, arma::field<arma::mat>&& W_B) {
  
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
  
  arma::mat At(S, S);
  arma::field<arma::mat> Bt(C);
  for (arma::uword c = 0; c < C; ++c) {
    Bt(c) = arma::mat(S, M(c));
  }
  arma::vec alpha(S);
  arma::vec alpha_new(S);
  arma::field<arma::vec> y_prob(C);
  arma::vec joint(arma::accu(M));
  for (arma::uword c = 0; c < C; ++c) {
    y_prob(c) = arma::vec(M(c));
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
      if (arma::all(obs(i).col(t - 1) < M)) {
        // previous observation is not missing, proceed normally
        state_prob(i).slice(t) = A.slice(t).each_col() % alpha;
        alpha = A.slice(t).t() * alpha;
        alpha_new.ones();
        for (arma::uword c = 0; c < C; ++c) {
          obs_prob(c, i).slice(t) = B(c).slice(t).cols(0, M(c) - 1).each_col() % alpha;
          if (obs(i)(c, t) < M(c)) {
            alpha_new %= obs_prob(c, i).slice(t).col(obs(i)(c, t));
          }
        }
      } else {
        // previous observation is missing, need to marginalize over it
        for (arma::uword c = 0; c < C; ++c) {
          if (obs(i)(c, t - 1) < M(c)) {
            y_prob(c).zeros();
            y_prob(obs(i)(c, t - 1)) = 1;
          } else {
            y_prob(c) = arma::sum(obs_prob(c, i).slice(t - 1)).t();
          }
        }
        joint = joint_probability(y_prob);
        for (arma::uword c = 0; c < C; ++c) {
          obs_prob(c, i).slice(t).zeros();
        }
        state_prob(i).slice(t).zeros();
        for (arma::uword j = 0; j < joint.n_elem; ++j) {
          for (arma::uword s = 0; s < S; ++s) {
            At.row(s) = softmax(gamma_A.slice(s) * W_A(j, i).col(t)).t();
          }
          state_prob(i).slice(t) += At.each_col() % alpha * joint(j);
          alpha_new = At.t() * alpha;
          for (arma::uword c = 0; c < C; ++c) {
            for (arma::uword s = 0; s < S; ++s) {
              Bt(c).row(s) = softmax(
                gamma_B(c).slice(s) * W_B(j, c, i).col(t)
              ).t();
            }
            obs_prob(c, i).slice(t) += Bt(c).each_col() % alpha_new * joint(j);
          }
        }
      }
      alpha_new.ones();
      for (arma::uword c = 0; c < C; ++c) {
        if (obs(i)(c, t) < M(c)) {
          alpha_new %= obs_prob(c, i).slice(t).col(obs(i)(c, t));
        } else {
          alpha_new %= arma::sum(obs_prob(c, i).slice(t), 1);
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

Rcpp::List fanhmm::simulate(
    arma::field<arma::mat>&& W_A, arma::field<arma::mat>&& W_B) {
  arma::field<arma::umat> y(N);
  arma::field<arma::uvec> z(N);
  arma::uvec seqS = arma::linspace<arma::uvec>(0, S - 1, S);
  arma::field<arma::uvec> seqM(C);
  arma::field<arma::vec> B(C);
  for (arma::uword c = 0; c < C; ++c) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
    B(c) = arma::vec(M(c));
  }
  arma::vec pi(S);
  arma::vec A(S);
  for (arma::uword i = 0; i < N; ++i) {
    y(i) = arma::umat(C, Ti(i));
    z(i) = arma::uvec(Ti(i));
    if (!icpt_only_pi || i == 0) {
      pi = softmax(gamma_pi * X_pi.col(i));
    }
    z(i)(0) = sample(seqS, pi.t());
    if (fixed_0) {
      for (arma::uword c = 0; c < C; ++c) {
        B(c) = softmax(
          gamma_B(c).slice(z(i)(0)) * X_B(c, i).col(0)
        );
        y(i)(c, 0) = sample(seqM(c), B(c).t());
      }
    } else {
      for (arma::uword c = 0; c < C; ++c) {
        B(c).zeros();
        for (arma::uword i = 0; prior_y.n_elem; ++i) {
          for (arma::uword m = 0; m < M(c); ++m) { // time
              B(c) += softmax(
                gamma_B(c).slice(z(i)(0)) * W_X_B(c, i).col(m)
              ) * prior_y(i);
          }
        }
        y(i)(c, 0) = sample(seqM(c), B(c).t());
      }
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      arma::uword idx = get_W_idx(y(i).col(t - 1), M);
      A = softmax(
        gamma_A.slice(z(i)(t - 1)) * W_A(idx, i).col(t)
      );
      z(i)(t) = sample(seqS, A.t());
      for (arma::uword c = 0; c < C; ++c) {
        B(c) = softmax(
          gamma_B(c).slice(z(i)(t)) * W_B(idx, c, i).col(t)
        );
        y(i)(c, t) = sample(seqM(c), B(c).t());
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y),
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
void fanhmm::gradient_B_t1(
    arma::mat& grad,
    arma::vec& tmpvec,
    const arma::mat& log_beta,
    const double ll,
    const arma::uword i,
    const arma::uword s,
    const arma::uword c) {
  
  arma::uword C = M.n_elem;
  arma::rowvec Brow(M(c));
  arma::uword idx = obs(i)(c, 0);
  if (!fixed_0) {
    double tmp = log_pi(s) + log_beta(s, 0) - ll;
    for (arma::uword j = 0; j < prior_y.n_elem; ++j) {
      Brow = B1(c).slice(j).row(s).cols(0, M(c) - 1);
      double brow = Brow(idx);
      tmpvec = -Brow.t() * brow;
      tmpvec(idx) += brow;
      double logpy = 0;
      for (arma::uword cc = 0; cc < C; ++cc) {
        if (cc != c) {
          logpy += log_B1(cc)(s, obs(i)(cc, 0), j);
        }
      }
      grad += exp(tmp + logpy) * tmpvec * W_X_B(j, c, i).t() * prior_y(j);
    }
  } else {
    Brow = B(c).slice(0).row(s).cols(0, M(c) - 1);
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
}