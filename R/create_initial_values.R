
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
    create_gamma_B_raw_nhmm(x[(n[i] + 1):(n[i + 1])], S, M[i], K)
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
  n <- sum((M - 1) * K * S)
  lapply(seq_len(D), function(i) {
    create_gamma_multichannel_B_raw_nhmm(x[(i - 1) * n + 1:n], S, M, K)
  })
}
create_gamma_omega_raw_mnhmm <- function(x, D, K) {
  matrix(x, D - 1, K)
}

create_gamma_pi_inits <- function(x, S, K, init_sd = 0, D = 1) {
  if (D > 1) {
    if (is.null(x)) {
      create_gamma_pi_raw_mnhmm(rnorm((S - 1) * K * D, sd = init_sd), S, K, D)
    } else {
      stopifnot_(
        length(unlist(x)) == (S - 1) * K * D,
        paste0(
          "Number of initial values for {.val gamma_pi} is not equal to ",
          "(S - 1) * K * D = {(S - 1) * K * D}."
        )
      )
      create_gamma_pi_raw_mnhmm(x, S, K, D)
    }
  } else {
    if (is.null(x)) {
      create_gamma_pi_raw_nhmm(rnorm((S - 1) * K, sd = init_sd), S, K)
    } else {
      stopifnot_(
        length(x) == (S - 1) * K,
        paste0(
          "Number of initial values for {.val gamma_pi} is not equal to ",
          "(S - 1) * K = {(S - 1) * K}."
        )
      )
      create_gamma_pi_raw_nhmm(x, S, K)
    }
  }
}
create_gamma_A_inits <- function(x, S, K, init_sd = 0, D = 1) {
  if (D > 1) {
    if (is.null(x)) {
      create_gamma_A_raw_mnhmm(rnorm((S - 1) * K * S * D, sd = init_sd), S, K, D)
    } else {
      stopifnot_(
        length(unlist(x)) == (S - 1) * K * S * D,
        paste0(
          "Number of initial values for {.val gamma_A} is not equal to ",
          "(S - 1) * K * S * D = {(S - 1) * K * S * D}."
        )
      )
      create_gamma_A_raw_mnhmm(x, S, K, D)
    }
  } else {
    if (is.null(x)) {
      create_gamma_A_raw_nhmm(rnorm((S - 1) * K * S, sd = init_sd), S, K)
    } else {
      stopifnot_(
        length(x) == (S - 1) * K * S,
        paste0(
          "Number of initial values for {.val gamma_A} is not equal to ",
          "(S - 1) * K * S = {(S - 1) * K * S}."
        )
      )
      create_gamma_A_raw_nhmm(x, S, K)
    }
  }
}
create_gamma_B_inits <- function(x, S, M, K, init_sd = 0, D = 1) {
  if (length(M) > 1) {
    if (D > 1) {
      if (is.null(x)) {
        create_gamma_multichannel_B_raw_mnhmm(
          rnorm(sum((M - 1) * K * S) * D, sd = init_sd), S, M, K, D
        )
      } else {
        stopifnot_(
          length(unlist(x)) == sum((M - 1) * K * S) * D,
          paste0(
            "Number of initial values for {.val gamma_B} is not equal to ",
            "sum((M - 1) * K * S) * D = {sum((M - 1) * K * S) * D}."
          )
        )
        create_gamma_multichannel_B_raw_mnhmm(x, S, M, K, D)
      }
    } else {
      if (is.null(x)) {
        create_gamma_multichannel_B_raw_nhmm(
          rnorm(sum((M - 1) * K * S), sd = init_sd), S, M, K
        )
      } else {
        stopifnot_(
          length(x) == sum((M - 1) * K * S),
          paste0(
            "Number of initial values for {.val gamma_B} is not equal to ",
            "sum((M - 1) * K * S) = {sum((M - 1) * K * S)}."
          )
        )
        create_gamma_multichannel_B_raw_nhmm(x, S, M, K)
      }
    }
  } else {
    if (D > 1) {
      if (is.null(x)) {
        create_gamma_B_raw_mnhmm(
          rnorm((M - 1) * K * S * D, sd = init_sd), S, M, K, D
        )
      } else {
        stopifnot_(
          length(unlist(x)) == (M - 1) * K * S * D,
          paste0(
            "Number of initial values for {.val gamma_B} is not equal to ",
            "(M - 1) * K * S * D = {(M - 1) * K * S * D}."
          )
        )
        create_gamma_B_raw_mnhmm(x, S, M, K, D)
      }
    } else {
      if (is.null(x)) {
        create_gamma_B_raw_nhmm(rnorm((M - 1) * K * S, sd = init_sd), S, M, K)
      } else {
        stopifnot_(
          length(x) == (M - 1) * K * S,
          paste0(
            "Number of initial values for {.val gamma_B} is not equal to ",
            "(M - 1) * K * S = {(M - 1) * K * S}."
          )
        )
        create_gamma_B_raw_nhmm(x, S, M, K)
      }
    }
  }
}
create_gamma_omega_inits <- function(x, D, K, init_sd = 0) {
  
  if (is.null(x)) {
    create_gamma_omega_raw_mnhmm(rnorm((D - 1) * K, sd = init_sd), D, K)
  } else {
    stopifnot_(
      length(x) == (D - 1) * K,
      paste0(
        "Number of initial values for {.val gamma_omega} is not equal to ",
        "(D - 1) * K = {(D - 1) * K}."
      )
    )
    create_gamma_omega_raw_mnhmm(x, D, K)
  }
}
#' Convert Initial Values for Inverse Softmax Scale
#' @noRd
create_inits_vector <- function(x, n, K, sd = 0, D = 1) {
  cbind(
    inv_softmax(x)[-1], 
    matrix(rnorm((n - 1) * (K - 1), sd = sd), n - 1, K - 1)
  )
}
create_inits_matrix <- function(x, n, m, K, sd = 0) {
  z <- array(0, c(m - 1, K, n))
  for (i in seq_len(n)) {
    z[, , i] <- t(create_inits_vector(x[i, ], m, K, sd))
  }
  z
}

create_initial_values <- function(inits, S, M, init_sd, K_i, K_s, K_o, K_d = 0, 
                                  D = 1) {
  
  if(!is.null(inits$initial_probs)) {
    if (D > 1) {
      gamma_pi_raw <- lapply(
        seq_len(D), function(i) {
          create_inits_vector(inits$initial_probs[[i]], S, K_i, init_sd)
        }
      )
    } else {
      gamma_pi_raw <- create_inits_vector(inits$initial_probs, S, K_i, init_sd)
    }
  } else {
    gamma_pi_raw <- create_gamma_pi_inits(
      inits$gamma_pi, S, K_i, init_sd, D
    )
  }
  
  if(!is.null(inits$transition_probs)) {
    if (D > 1) {
      gamma_A_raw <- lapply(
        seq_len(D), function(i) {
          create_inits_matrix(inits$transition_probs[[i]], S, S, K_s, init_sd)
        }
      )
    } else {
      gamma_A_raw <- create_inits_matrix(
        inits$transition_probs, S, S, K_s, init_sd
      )
    }
  } else {
    gamma_A_raw <- create_gamma_A_inits(inits$gamma_A, S, K_s, init_sd, D)
  }
  
  if(!is.null(inits$emission_probs)) {
    if (D > 1) {
      if (length(M) > 1) {
        gamma_B_raw <-  lapply(
          seq_len(D), function(i) {
            lapply(seq_len(length(M)), function(j) {
              create_inits_matrix(
                inits$emission_probs[[i]][[j]], S, M[j], K_o, init_sd)
            })
          })
      } else {
        gamma_B_raw <- lapply(
          seq_len(D), function(i) {
            create_inits_matrix(inits$emission_probs[[i]], S, M, K_o, init_sd)
          }
        )
      }
    } else {
      if (length(M) > 1) {
        gamma_B_raw <- lapply(seq_len(length(M)), function(j) {
          create_inits_matrix(
            inits$emission_probs[[j]], S, M[j], K_o, init_sd)
        })
      } else {
        gamma_B_raw <- create_inits_matrix(
          inits$emission_probs, S, M, K_o, init_sd
        )
      }
    }
  } else {
    gamma_B_raw <- create_gamma_B_inits(inits$gamma_B, S, M, K_o, init_sd, D)
  }
  out <- list(
    gamma_pi_raw = gamma_pi_raw,
    gamma_A_raw = gamma_A_raw,
    gamma_B_raw = gamma_B_raw
  )
  if (D > 1) {
    if(!is.null(inits$cluster_probs)) {
      gamma_omega_raw <- create_inits_vector(
        inits$cluster_probs, D, K_d, init_sd
      )
    } else {
      gamma_omega_raw <- create_gamma_omega_inits(
        inits$gamma_omega, D, K_d, init_sd
      )
    }
    out$gamma_omega_raw <- gamma_omega_raw
  }
  out
}
