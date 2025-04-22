#include "config.h"
#include "fanhmm.h"
#include "softmax.h"
#include "list_to_field.h"
#include "joint_probability.h"
#include "sample.h"

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
    W_X_B(fixed_0 ? arma::field<arma::vec>() : veclist_to_3d_field(W_X_B)) {
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
          Btmp.rows(0, M(c) - 1).zeros();
          for (arma::uword i = 0; prior_y.n_elem; ++i) {
            for (arma::uword m = 0; m < M(c); ++m) { // time
              for (arma::uword s = 0; s < S; ++s) { // from states
                Btmp.col(s).rows(0, M(c) - 1) += softmax(
                  gamma_B(c).slice(s) * W_X_B(c, i).col(m)
                ) * prior_y(i);
              }
            }
          }
          B(c).slice(0) = Btmp.t();
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
  for (arma::uword c = 0; c < C; ++c) {
    seqM(c) = arma::linspace<arma::uvec>(0, M(c) - 1, M(c));
  }
  arma::vec pi(S);
  arma::vec A(S);
  arma::field<arma::vec> B(C);
  for (arma::uword i = 0; i < N; ++i) {
    y(i) = arma::umat(C, Ti(i));
    z(i) = arma::uvec(Ti(i));
    if (!icpt_only_pi || i == 0) {
      pi = softmax(gamma_pi * X_pi.col(i));
    }
    z(i)(0) = sample(seqS, pi);
    if (fixed_0) {
      for (arma::uword c = 0; c < C; ++c) {
        B(c) = softmax(
          gamma_B(c).slice(z(i)(0)) * X_B(c, i).col(0)
        );
        y(i)(c, 0) = sample(seqM(c), B(c));
      }
    } else {
      for (arma::uword c = 0; c < C; ++c) {
        arma::mat Btmp(M(c) + 1, S, arma::fill::ones);
        Btmp.rows(0, M(c) - 1).zeros();
        for (arma::uword i = 0; prior_y.n_elem; ++i) {
          for (arma::uword m = 0; m < M(c); ++m) { // time
            for (arma::uword s = 0; s < S; ++s) { // from states
              Btmp.col(s).rows(0, M(c) - 1) += softmax(
                gamma_B(c).slice(s) * W_X_B(c, i).col(m)
              ) * prior_y(i);
            }
          }
        }
        B(c) = Btmp.t();
        y(i)(c, 0) = sample(seqM(c), B(c));
      }
    }
    for (arma::uword t = 1; t < Ti(i); ++t) {
      for (arma::uword c = 0; c < C; ++c) {
        A = softmax(
          gamma_A.slice(z(i)(t - 1)) * W_A(y(i)(c, t - 1), i).col(t)
        );
        z(i)(t) = sample(seqS, A);
        B(c) = softmax(
          gamma_B(c).slice(z(i)(t)) * W_B(y(i)(c, t - 1)).col(t)
        );
        y(i)(c, t) = sample(seqM(c), B(c));
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("observations") = Rcpp::wrap(y),
    Rcpp::Named("states") = Rcpp::wrap(z)
  );
}
