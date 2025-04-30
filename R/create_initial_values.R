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
  K_pi <- K(model$X_pi)
  K_A <- K(model$X_A)
  K_B <- K(model$X_B)
  if (D > 1) {
    K_omega <- K(model$X_omega)
  } else {
    K_omega <- 0
    D <- 1
  }
  create_initial_values_(inits, init_sd, S, M, K_pi, K_A, K_B, D, K_omega)
}
create_initial_values_ <- function(inits, init_sd, S, M, K_pi, K_A, K_B, D = 1, 
                                   K_omega = 0) {
  if(!is.null(inits$initial_probs)) {
    if (D > 1) {
      stopifnot_(
        is.list(inits$initial_probs) && length(inits$initial_probs) == D,
        "Initial values for initial states of MNHMM must be a list of length D,
        with each element containing the initial state probability vector of 
        length S, where D is the nubmer of mixtures and S is the number of 
        states."
      )
      eta_pi <- lapply(
        seq_len(D), \(i) {
          pi <- inits$initial_probs[[i]]
          stopifnot_(
            length(pi) == S,
            "Initial values for initial states of MNHMM must be a list of length D,
            with each element containing the initial state probability vector of 
            length S, where D is the nubmer of mixtures and S is the number of 
            states."
          )
          create_inits_vector(pi, S, K_pi, init_sd)
        }
      )
    } else {
      pi <- inits$initial_probs
      stopifnot_(
        length(pi) == S,
        "Initial values for initial state probabilities must be a vector of 
        length S, matching the number of states."
      )
      eta_pi <- create_inits_vector(pi, S, K_pi, init_sd)
    }
  } else {
    eta_pi <- create_eta_pi_inits(
      inits$eta_pi, S, K_pi, init_sd, D
    )
  }
  
  if(!is.null(inits$transition_probs)) {
    if (D > 1) {
      stopifnot_(
        is.list(inits$transition_probs) && length(inits$transition_probs) == D,
        "Initial values for transitions of MNHMM must be a list of length D, 
        with each element containing S x S transition matrices, where 
        D is the number of mixtures and S is the number of states."
      )
      eta_A <- lapply(
        seq_len(D), \(i) {
          A <- inits$transition_probs[[i]]
          stopifnot_(
            nrow(A) == S && ncol(A) == S,
            "Initial transition probability matrices must have dimensions 
              (number of states) x (number of states)."
          )
          create_inits_matrix(A, S, S, K_A, init_sd)
        }
      )
    } else {
      A <- inits$transition_probs
      stopifnot_(
        nrow(A) == S && ncol(A) == S,
        "Initial transition probability matrix must have dimensions 
        (number of states) x (number of states)."
      )
      eta_A <- create_inits_matrix(A, S, S, K_A, init_sd)
    }
  } else {
    eta_A <- create_eta_A_inits(inits$eta_A, S, K_A, init_sd, D)
  }
  
  if(!is.null(inits$emission_probs)) {
    C <- length(M)
    if (D > 1) {
      stopifnot_(
        is.list(inits$emission_probs) && length(inits$emission_probs) == D,
        "Initial values for emissions of MNHMM must be a list of length D, 
        where each element contains C matrices, one for each response."
      )
      eta_B <- lapply(
        seq_len(D), \(i) {
          if (C == 1 && is.matrix(inits$emission_probs[[i]])) {
            inits$emission_probs[[i]] <- list(inits$emission_probs[[i]])
          }
          stopifnot_(
            is.list(inits$emission_probs[[i]]) && length(inits$emission_probs[[i]]) == C,
            "Initial values for emissions of MNHMM must be a list of length D, 
            where each element contains C matrices, one for each response."
          )
          lapply(seq_along(M), \(j) {
            B <- inits$emission_probs[[i]][[j]]
            stopifnot_(
              nrow(B) == S && ncol(B) == M[j],
              "Initial emission probability matrices must have dimensions 
              (number of states) x (number of symbols)."
            )
            create_inits_matrix(B, S, M[j], K_B[j], init_sd)
          })
        })
    } else {
      if (C == 1 && is.matrix(inits$emission_probs)) {
        inits$emission_probs <- list(inits$emission_probs)
      }
      stopifnot_(
        is.list(inits$emission_probs) && length(inits$emission_probs) == C,
        "Initial values for emissions of NHMM must be a list of length C, 
        containing the emission probabilties for each response."
      )
      eta_B <- lapply(seq_along(M), \(j) {
        B <- inits$emission_probs[[j]]
        stopifnot_(
          nrow(B) == S && ncol(B) == M[j],
          "Initial emission probability matrices must have dimensions 
              (number of states) x (number of symbols)."
        )
        create_inits_matrix(B, S, M[j], K_B[j], init_sd)
      })
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
      stopifnot_(
        length(inits$cluster_probs) == S,
        "Initial cluster probabilities must be a vector with length D, 
        matching the number of mixtures"
      )
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
  if (sd > 0) {
    x <- stats::rnorm((S - 1) * K, x, sd)
  }
  matrix(x, S - 1, K)
}
create_eta_A_nhmm <- function(x, S, K, sd = 0) {
  if (sd > 0) {
    x <- stats::rnorm((S - 1) * K * S, x, sd)
  }
  array(x, c(S - 1, K, S))
}
create_eta_B_nhmm <- function(x, S, M, K, sd = 0) {
  n <- c(0, cumsum((M - 1) * K * S))
  if (sd > 0) {
    x <- stats::rnorm(n[length(n)], x, sd)
  }
  lapply(seq_along(M), \(i) {
    array(x[(n[i] + 1):(n[i + 1])], c(M[i] - 1, K[i], S))
  })
}
create_eta_pi_mnhmm <- function(x, S, K, D, sd = 0) {
  n <- (S - 1) * K
  lapply(seq_len(D), \(i) {
    create_eta_pi_nhmm(x[(i - 1) * n + 1:n], S, K, sd)
  })
}
create_eta_A_mnhmm <- function(x, S, K, D, sd = 0) {
  n <- (S - 1) * K * S
  lapply(seq_len(D), \(i) {
    create_eta_A_nhmm(x[(i - 1) * n + 1:n], S, K, sd)
  })
}

create_eta_B_mnhmm <- function(x, S, M, K, D, sd = 0) {
  n <- sum((M - 1) * K * S)
  lapply(seq_len(D), \(i) {
    create_eta_B_nhmm(x[(i - 1) * n + 1:n], S, M, K, sd)
  })
}
create_eta_omega_mnhmm <- function(x, D, K, sd = 0) {
  matrix(stats::rnorm((D - 1) * K, x, sd), D - 1, K)
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
  
  if (D > 1) {
    if (is.null(x)) {
      create_eta_B_mnhmm(
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
      create_eta_B_mnhmm(unlist(x), S, M, K, D, init_sd)
    }
  } else {
    if (is.null(x)) {
      create_eta_B_nhmm(
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
      create_eta_B_nhmm(unlist(x), S, M, K, init_sd)
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
    matrix(stats::rnorm((n - 1) * (K - 1), sd = sd), n - 1, K - 1)
  )
}
create_inits_matrix <- function(x, n, m, K, sd = 0) {
  z <- array(0, c(m - 1, K, n))
  for (i in seq_len(n)) {
    z[, , i] <- create_inits_vector(x[i, ], m, K, sd)
  }
  z
}
