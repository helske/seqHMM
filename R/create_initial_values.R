#' Generate Initial values for NHMMs
#' 
#' `create_initial_values` generates initial values for the `eta` parameters 
#' (regression coefficients) for the NHMM and MNHMMs.
#' 
#' @param inits List (possibly empty).
#' @noRd
create_initial_values <- function(inits, model, init_sd) {
  S <- model$n_states
  M <- model$n_symbols
  D <- model$n_clusters
  K_pi <- nrow(model$X_pi)
  K_A <- nrow(model$X_A)
  K_B <- nrow(model$X_B)
  if (D > 1) {
    K_omega <- nrow(model$X_omega)
  } else {
    D <- 1
  }
  create_initial_values_(inits, init_sd, S, M, K_pi, K_A, K_B, D, K_omega)
}
create_initial_values_ <- function(inits, init_sd, S, M, K_pi, K_A, K_B, D = 1, 
                                   K_omega = 0) {
  if(!is.null(inits$initial_probs)) {
    if (D > 1) {
      eta_pi <- lapply(
        seq_len(D), function(i) {
          create_inits_vector(inits$initial_probs[[i]], S, K_pi, init_sd)
        }
      )
    } else {
      eta_pi <- create_inits_vector(inits$initial_probs, S, K_pi, init_sd)
    }
  } else {
    eta_pi <- create_eta_pi_inits(
      inits$eta_pi, S, K_pi, init_sd, D
    )
  }
  
  if(!is.null(inits$transition_probs)) {
    if (D > 1) {
      eta_A <- lapply(
        seq_len(D), function(i) {
          create_inits_matrix(inits$transition_probs[[i]], S, S, K_A, init_sd)
        }
      )
    } else {
      eta_A <- create_inits_matrix(
        inits$transition_probs, S, S, K_A, init_sd
      )
    }
  } else {
    eta_A <- create_eta_A_inits(inits$eta_A, S, K_A, init_sd, D)
  }
  
  if(!is.null(inits$emission_probs)) {
    if (D > 1) {
      if (length(M) > 1) {
        eta_B <- lapply(
          seq_len(D), function(i) {
            lapply(seq_len(length(M)), function(j) {
              create_inits_matrix(
                inits$emission_probs[[i]][[j]], S, M[j], K_B, init_sd)
            })
          })
      } else {
        eta_B <- lapply(
          seq_len(D), function(i) {
            create_inits_matrix(inits$emission_probs[[i]], S, M, K_B, init_sd)
          }
        )
      }
    } else {
      if (length(M) > 1) {
        eta_B <- lapply(seq_len(length(M)), function(j) {
          create_inits_matrix(
            inits$emission_probs[[j]], S, M[j], K_B, init_sd)
        })
      } else {
        eta_B <- create_inits_matrix(
          inits$emission_probs, S, M, K_B, init_sd
        )
      }
    }
  } else {
    eta_B <- create_eta_B_inits(inits$eta_B, S, M, K_B, init_sd, D)
  }
  out <- list(
    eta_pi = eta_pi,
    eta_A = eta_A,
    eta_B = eta_B
  )
  if (D > 1) {
    if(!is.null(inits$cluster_probs)) {
      eta_omega <- create_inits_vector(
        inits$cluster_probs, D, K_omega, init_sd
      )
    } else {
      eta_omega <- create_eta_omega_inits(
        inits$eta_omega, D, K_omega, init_sd
      )
    }
    out$eta_omega <- eta_omega
  }
  out
}

create_eta_pi_nhmm <- function(x, S, K, sd = 0) {
  matrix(rnorm((S - 1) * K, x, sd), S - 1, K)
}
create_eta_A_nhmm <- function(x, S, K, sd = 0) {
  array(rnorm((S - 1) * K * S, x, sd), c(S - 1, K, S))
}
create_eta_B_nhmm <- function(x, S, M, K, sd = 0) {
  array(rnorm((M - 1) * K * S, x, sd), c(M - 1, K, S))
  
}
create_eta_multichannel_B_nhmm <- function(x, S, M, K, sd = 0) {
  n <- c(0, cumsum((M - 1) * K * S))
  lapply(seq_len(length(M)), function(i) {
    create_eta_B_nhmm(x[(n[i] + 1):(n[i + 1])], S, M[i], K, sd)
  })
}
create_eta_pi_mnhmm <- function(x, S, K, D, sd = 0) {
  n <- (S - 1) * K
  lapply(seq_len(D), function(i) {
    create_eta_pi_nhmm(x[(i - 1) * n + 1:n], S, K, sd)
  })
}
create_eta_A_mnhmm <- function(x, S, K, D, sd = 0) {
  n <- (S - 1) * K * S
  lapply(seq_len(D), function(i) {
    create_eta_A_nhmm(x[(i - 1) * n + 1:n], S, K, sd)
  })
}
create_eta_B_mnhmm <- function(x, S, M, K, D, sd = 0) {
  n <- (M - 1) * K * S
  lapply(seq_len(D), function(i) {
    create_eta_B_nhmm(x[(i - 1) * n + 1:n], S, M, K, sd)
  })
}
create_eta_multichannel_B_mnhmm <- function(x, S, M, K, D, sd = 0) {
  n <- sum((M - 1) * K * S)
  lapply(seq_len(D), function(i) {
    create_eta_multichannel_B_nhmm(x[(i - 1) * n + 1:n], S, M, K, sd)
  })
}
create_eta_omega_mnhmm <- function(x, D, K, sd = 0) {
  matrix(rnorm((D - 1) * K, x, sd), D - 1, K)
}

create_eta_pi_inits <- function(x, S, K, init_sd = 0, D = 1) {
  if (D > 1) {
    if (is.null(x)) {
      create_eta_pi_mnhmm(numeric((S - 1) * K * D), S, K, D, init_sd)
    } else {
      stopifnot_(
        length(unlist(x)) == (S - 1) * K * D,
        paste0(
          "Number of initial values for {.val eta_pi} is not equal to ",
          "(S - 1) * K * D = {(S - 1) * K * D}."
        )
      )
      create_eta_pi_mnhmm(unlist(x), S, K, D, init_sd)
    }
  } else {
    if (is.null(x)) {
      create_eta_pi_nhmm(numeric((S - 1) * K), S, K, init_sd)
    } else {
      stopifnot_(
        length(x) == (S - 1) * K,
        paste0(
          "Number of initial values for {.val eta_pi} is not equal to ",
          "(S - 1) * K = {(S - 1) * K}."
        )
      )
      create_eta_pi_nhmm(x, S, K, init_sd)
    }
  }
}
create_eta_A_inits <- function(x, S, K, init_sd = 0, D = 1) {
  if (D > 1) {
    if (is.null(x)) {
      create_eta_A_mnhmm(numeric((S - 1) * K * S * D), S, K, D, init_sd)
    } else {
      stopifnot_(
        length(unlist(x)) == (S - 1) * K * S * D,
        paste0(
          "Number of initial values for {.val eta_A} is not equal to ",
          "(S - 1) * K * S * D = {(S - 1) * K * S * D}."
        )
      )
      create_eta_A_mnhmm(unlist(x), S, K, D, init_sd)
    }
  } else {
    if (is.null(x)) {
      create_eta_A_nhmm(numeric((S - 1) * K * S), S, K, init_sd)
    } else {
      stopifnot_(
        length(x) == (S - 1) * K * S,
        paste0(
          "Number of initial values for {.val eta_A} is not equal to ",
          "(S - 1) * K * S = {(S - 1) * K * S}."
        )
      )
      create_eta_A_nhmm(x, S, K, init_sd)
    }
  }
}
create_eta_B_inits <- function(x, S, M, K, init_sd = 0, D = 1) {
  if (length(M) > 1) {
    if (D > 1) {
      if (is.null(x)) {
        create_eta_multichannel_B_mnhmm(
          numeric(sum((M - 1) * K * S) * D), 
          S, M, K, D, init_sd
        )
      } else {
        stopifnot_(
          length(unlist(x)) == sum((M - 1) * K * S) * D,
          paste0(
            "Number of initial values for {.val eta_B} is not equal to ",
            "sum((M - 1) * K * S) * D = {sum((M - 1) * K * S) * D}."
          )
        )
        create_eta_multichannel_B_mnhmm(unlist(x), S, M, K, D, init_sd)
      }
    } else {
      if (is.null(x)) {
        create_eta_multichannel_B_nhmm(
          numeric(sum((M - 1) * K * S)), S, M, K, init_sd
        )
      } else {
        stopifnot_(
          length(unlist(x)) == sum((M - 1) * K * S),
          paste0(
            "Number of initial values for {.val eta_B} is not equal to ",
            "sum((M - 1) * K * S) = {sum((M - 1) * K * S)}."
          )
        )
        create_eta_multichannel_B_nhmm(unlist(x), S, M, K, init_sd)
      }
    }
  } else {
    if (D > 1) {
      if (is.null(x)) {
        create_eta_B_mnhmm(numeric((M - 1) * K * S * D), S, M, K, D, init_sd)
      } else {
        stopifnot_(
          length(unlist(x)) == (M - 1) * K * S * D,
          paste0(
            "Number of initial values for {.val eta_B} is not equal to ",
            "(M - 1) * K * S * D = {(M - 1) * K * S * D}."
          )
        )
        create_eta_B_mnhmm(unlist(x), S, M, K, D, init_sd)
      }
    } else {
      if (is.null(x)) {
        create_eta_B_nhmm(numeric((M - 1) * K * S), S, M, K, init_sd)
      } else {
        stopifnot_(
          length(x) == (M - 1) * K * S,
          paste0(
            "Number of initial values for {.val eta_B} is not equal to ",
            "(M - 1) * K * S = {(M - 1) * K * S}."
          )
        )
        create_eta_B_nhmm(x, S, M, K, init_sd)
      }
    }
  }
}
create_eta_omega_inits <- function(x, D, K, init_sd = 0) {
  
  if (is.null(x)) {
    create_eta_omega_mnhmm(numeric((D - 1) * K), D, K, init_sd)
  } else {
    stopifnot_(
      length(x) == (D - 1) * K,
      paste0(
        "Number of initial values for {.val eta_omega} is not equal to ",
        "(D - 1) * K = {(D - 1) * K}."
      )
    )
    create_eta_omega_mnhmm(x, D, K, init_sd)
  }
}
#' Convert Initial Values for Inverse Softmax Scale
#' @noRd
create_inits_vector <- function(x, n, K, sd = 0, D = 1) {
  cbind(
    p_to_eta(x), # intercepts
    matrix(rnorm((n - 1) * (K - 1), sd = sd), n - 1, K - 1)
  )
}
create_inits_matrix <- function(x, n, m, K, sd = 0) {
  z <- array(0, c(m - 1, K, n))
  for (i in seq_len(n)) {
    z[, , i] <- create_inits_vector(x[i, ], m, K, sd)
  }
  z
}

create_rho_A_inits <- function(x, S, M, L, init_sd = 0) {
  if (is.null(x$rho_A)) {
    create_rho_A(numeric((S - 1) * L * (M - 1) * S), S, M, L, init_sd)
  } else {
    stopifnot_(
      length(x$rho_A) == (S - 1) * L * (M - 1) * S,
      paste0(
        "Number of initial values for {.val rho_A} is not equal to ",
        "(S - 1) * L * (M - 1) * S = {(S - 1) * L * (M - 1) * S}."
      )
    )
    create_rho_A(x$rho_A, S, M, L, init_sd)
  }
}

create_rho_A <- function(x, S, M, L, init_sd = 0) {
  n <- (S - 1) * L * (M - 1)
  lapply(seq_len(S), function(i) {
    array(rnorm(n, x[(i - 1) * n + 1:n], init_sd), c(S - 1, L, M - 1))
  })
}

create_rho_B_inits <- function(x, S, M, L, init_sd = 0) {
  if (is.null(x$rho_B)) {
    create_rho_B(numeric((M - 1) * L * (M - 1) * S), S, M, L, init_sd)
  } else {
    stopifnot_(
      length(x$rho_B) == (M - 1) * L * (M - 1) * S,
      paste0(
        "Number of initial values for {.val rho_B} is not equal to ",
        "(M - 1) * L * (M - 1) * S = {(M - 1) * L * (M - 1) * S}."
      )
    )
    create_rho_B(x$rho_B, S, M, L, init_sd)
  }
}

create_rho_B <- function(x, S, M, L, init_sd = 0) {
  n <- (M - 1) * L * (M - 1)
  lapply(seq_len(S), function(i) {
    array(rnorm(n, x[(i - 1) * n + 1:n], init_sd), c(M - 1, L, M - 1))
  })
}

