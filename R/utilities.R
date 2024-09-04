#' Filter data.frame Rows with Time not Exceeding Sequence Length
#' @noRd
non_void_rows <- function(d, time, id, sequence_lengths) {
  d[d[[time]] <= sequence_lengths[match(d[[id]], unique(d[[id]]))], ]
}
#' Check If an Object Is a List of Lists
#' 
#' @noRd
is_list_of_lists <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  } else {
    if (all(unlist(lapply(x, is.list)))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
#' Regularized Inverse Softmax Function
#' 
#' @noRd
inv_softmax <- function(x) {
  x <- pmin(pmax(x, 0.001), 0.999)
  x <- x / sum(x)
  log(x) - log(x[1])
}
#' Bind list of 3D Arrays to a 4D Array
#' @noRd
bind_3D_arrays <- function(x) {
  y <- array(0, c(length(x), dim(x[[1]])))
  for (i in seq_along(x)) {
    y[i,,,] <- x[[i]]
  }
  y
}
#' Bind list of 2D Arrays to a 3D Array
#' @noRd
bind_2D_arrays <- function(x) {
  y <- array(0, c(length(x), dim(x[[1]])))
  for (i in seq_along(x)) {
    y[i,,] <- x[[i]]
  }
  y
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
    z <- array(rnorm(n * (m - 1) * K, sd = sd), c(n, m - 1, K))
  } else {
    z <- array(0, c(n, m - 1, K))
    for (i in seq_len(n)) {
      z[i, ,] <- create_inits_vector(x[i, ], m, K, sd)
    }
  }
  z
}
create_inits_multichannel <- function(x, n, m, K, sd = 0) {
  if (is.null(x)) {
    rnorm(n * (sum(m) - length(m)) * K, sd = sd)
  } else {
    unlist(lapply(seq_along(m), function(i) {
      aperm(create_inits_matrix(x[[i]], n, m[i], K, sd), c(2, 1, 3))
    }))
  }
}

create_initial_values <- function(
    inits, S, M, init_sd, K_i, K_s, K_o, K_d = 0, D = 0) {
  if (D > 0) {
    list(
      beta_i_raw = if (is.null(inits$beta_i_raw)) {
        bind_2D_arrays(lapply(seq_len(D), function(i) {
          create_inits_vector(inits$initial_probs[[i]], S, K_i, init_sd)
        }))
      } else inits$beta_i_raw,
      beta_s_raw = if (is.null(inits$beta_i_raw)) {
        bind_3D_arrays(lapply(seq_len(D), function(i) {
          create_inits_matrix(inits$transition_probs[[i]], S, S, K_s, init_sd)
        }))
      } else inits$beta_s_raw,
      beta_o_raw = if (is.null(inits$beta_o_raw)) {
        if (length(M) > 1) {
          do.call("rbind", lapply(seq_len(D), function(i) {
            create_inits_multichannel(inits$emission_probs[[i]], S, M, K_o, init_sd)
          }))
        } else {
          bind_3D_arrays(lapply(seq_len(D), function(i) {
            create_inits_matrix(inits$emission_probs[[i]], S, M, K_o, init_sd)
          }))
        }
      } else inits$beta_o_raw,
      theta_raw = if (is.null(inits$theta_raw)) {
        create_inits_vector(inits$cluster_probs, D, K_d, init_sd)
      } else inits$theta_raw
    )
  } else {
    list(
      beta_i_raw = if (is.null(inits$beta_i_raw)) {
        create_inits_vector(inits$initial_probs, S, K_i, init_sd)
      } else inits$beta_i_raw,
      beta_s_raw = if (is.null(inits$beta_s_raw)) {
        create_inits_matrix(inits$transition_probs, S, S, K_s, init_sd) 
      } else inits$beta_s_raw,
      beta_o_raw =  if (is.null(inits$beta_o_raw)) {
        if (length(M) > 1) {
          create_inits_multichannel(inits$emission_probs, S, M, K_o, init_sd)
        } else {
          create_inits_matrix(inits$emission_probs, S, M, K_o, init_sd) 
        }
      } else {
        inits$beta_o_raw
      }
    )
  }
}
#' Stop Function Execution Unless Condition Is True
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}
#' Give a Informational Message
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_inform()].
#' @param ... See [cli::cli_inform()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}
#' Issue a Warning
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_warn()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
warning_ <- function (message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}
#' Issue an Error
#' 
#' Function copied from the `dynamite` package.
#' 
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#' 
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}
#' Checks Whether the Formula Contains Only the Intercept Term
#' 
#' @param f A formula object.
#' @noRd
intercept_only <- function(f) {
  identical(deparse(update(f, 0 ~ .)), "0 ~ 1")
}
#' Create obsArray for Various C++ functions
#' 
#' @noRd
create_obsArray <- function(model) {
  obsArray <- array(
    0L, 
    c(model$n_sequences, model$length_of_sequences, model$n_channels))
  if (model$n_channels == 1) {
    model$observations <- list(model$observations)
  }
  for (i in 1:model$n_channels) {
    obsArray[, , i] <- unlist(lapply(model$observations[[i]], as.integer)) - 1L
    obsArray[, , i][obsArray[, , i] > model$n_symbols[i]] <- model$n_symbols[i]
    stopifnot_(
      sum(obsArray[, , i] < model$n_symbols[i]) > 0,
      "One channel contains only missing values, model is degenerate."
    )
  }
  aperm(obsArray)
}
#' Create emissionArray for Various C++ functions
#' 
#' @noRd
create_emissionArray <- function(model) {
  emissionArray <- array(1, c(model$n_states, max(model$n_symbols) + 1, model$n_channels))
  if (model$n_channels == 1) {
    model$emission_probs <- list(model$emission_probs)
  }
  for (i in 1:model$n_channels) {
    emissionArray[, 1:model$n_symbols[i], i] <- model$emission_probs[[i]]
  }
  emissionArray
}

#' Remove rows corresponding to void symbols
#' @noRd
remove_voids <- function(model, x) {
  
  x_id <- x[[model$id_variable]]
  x_time <- x[[model$time_variable]]
  if (model$n_channels == 1) {
    time <- colnames(model$observations)
    id <- rownames(model$observations)
  } else {
    time <- colnames(model$observations[[1]])
    id <- rownames(model$observations[[1]])
  }
  do.call(
    "rbind", 
    lapply(model$n_sequences, function(i) {
      idx <- which(
        x_id == id[i] & x_time %in% time[seq_len(model$sequence_lengths[i])]
      )
      x[idx, ]
    })
  )
}
