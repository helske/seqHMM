
create_gamma_pi_raw_nhmm <- function(x, S, K) {
  matrix(x, S - 1, K)
}
create_gamma_A_raw_nhmm <- function(x, S, K) {
  array(x, c(S - 1, K, S))
}
create_gamma_B_raw_nhmm <- function(x, S, M, K) {
  array(x, c(M - 1, K, S))
}
create_gamma_multichannel_B_raw_nhmm <- function(x, S, M, K) {
  n <- c(0, cumsum((M - 1) * K * S))
  lapply(seq_len(length(M)), function(i) {
    create_gamma_B_raw_nhmm(x[(i - 1) * n[i] + 1:((M[i] - 1) * K * S)], S, M[i], K)
  })
}
create_gamma_pi_raw_mnhmm <- function(x, S, K, D) {
  n <- (S - 1) * K
  lapply(seq_len(D), function(i) {
    create_gamma_pi_raw_nhmm(x[(i - 1) * n + 1:n], S, K)
  })
}
create_gamma_A_raw_mnhmm <- function(x, S, K, D) {
  n <- (S - 1) * K * S
  lapply(seq_len(D), function(i) {
    create_gamma_A_raw_nhmm(x[(i - 1) * n + 1:n], S, K)
  })
}
create_gamma_B_raw_mnhmm <- function(x, S, M, K, D) {
  n <- (M - 1) * K * S
  lapply(seq_len(D), function(i) {
    create_gamma_B_raw_nhmm(x[(i - 1) * n + 1:n], S, M, K)
  })
}
create_gamma_multichannel_B_raw_mnhmm <- function(x, S, M, K, D) {
  n <- (sum(M) - 1) * K * S
  lapply(seq_len(D), function(i) {
    create_gamma_multichannel_B_raw_nhmm(x[(i - 1) * n + 1:n], S, M, K)
  })
}
create_gamma_omega_raw_mnhmm <- function(x, K, D) {
  matrix(x, D - 1, K)
}


#' Convert Initial Values for Inverse Softmax Scale
#' @noRd
create_inits_vector <- function(x, n, K, sd = 0) {
  if (is.null(x)) {
    matrix(rnorm((n - 1) * K, sd = sd), n - 1, K)
  } else {
    cbind(
      inv_softmax(x)[-1], 
      matrix(rnorm((n - 1) * (K - 1), sd = sd), n - 1, K - 1)
    )
  }
}
create_inits_matrix <- function(x, n, m, K, sd = 0) {
  if (is.null(x)) {
    z <- array(rnorm(n * (m - 1) * K, sd = sd), c(m - 1, K, n))
  } else {
    z <- array(0, c(m - 1, K, n))
    for (i in seq_len(n)) {
      z[, , i] <- t(create_inits_vector(x[i, ], m, K, sd))
    }
  }
  z
}

create_initial_values <- function(
    inits, S, M, init_sd, K_i, K_s, K_o, K_d = 0, D = 1) {
  if (D > 1) {
    list(
      gamma_pi_raw = if (is.null(inits$gamma_pi_raw)) {
        lapply(seq_len(D), function(i) {
          create_inits_vector(inits$initial_probs[[i]], S, K_i, init_sd)
        })
      } else inits$gamma_pi_raw,
      gamma_A_raw = if (is.null(inits$gamma_pi_raw)) {
        lapply(seq_len(D), function(i) {
          create_inits_matrix(inits$transition_probs[[i]], S, S, K_s, init_sd)
        })
      } else inits$gamma_A_raw,
      gamma_B_raw = if (is.null(inits$gamma_B_raw)) {
        if (length(M) > 1) {
          lapply(seq_len(D), function(i) {
            lapply(seq_len(length(M)), function(j) {
              create_inits_matrix(
                inits$emission_probs[[i]][[j]], S, M[j], K_o, init_sd)
            })
          })
        } else {
          lapply(seq_len(D), function(i) {
            create_inits_matrix(inits$emission_probs[[i]], S, M, K_o, init_sd)
          })
        }
      } else inits$gamma_B_raw,
      gamma_omega_raw = if (is.null(inits$gamma_omega_raw)) {
        create_inits_vector(inits$cluster_probs, D, K_d, init_sd)
      } else inits$gamma_omega_raw
    )
  } else {
    list(
      gamma_pi_raw = if (is.null(inits$gamma_pi_raw)) {
        create_inits_vector(inits$initial_probs, S, K_i, init_sd)
      } else inits$gamma_pi_raw,
      gamma_A_raw = if (is.null(inits$gamma_A_raw)) {
        create_inits_matrix(inits$transition_probs, S, S, K_s, init_sd) 
      } else inits$gamma_A_raw,
      gamma_B_raw = if (is.null(inits$gamma_B_raw)) {
        if (length(M) > 1) {
          lapply(seq_len(length(M)), function(j) {
            create_inits_matrix(
              inits$emission_probs[[j]], S, M[j], K_o, init_sd)
          })
        } else {
          create_inits_matrix(inits$emission_probs, S, M, K_o, init_sd) 
        }
      } else {
        inits$gamma_B_raw
      }
    )
  }
}
