#' #' Predict method for non-homogeneous hidden Markov models
#' #' 
#' #' This is essentially same as `get_probs` but with option to return samples.
#' #' 
#' #' @param object A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' #' @param newdata Optional data frame which is used for prediction.
#' #' @param ... Ignored.
#' #' @noRd
#' predict.nhmm <- function(
#'     object, newdata = NULL, dontchange_colnames = FALSE, ...) {
#'   
#'   if (!is.null(newdata)) {
#'     time <- object$time_variable
#'     id <- object$id_variable
#'     stopifnot_(
#'       is.data.frame(newdata), 
#'       "Argument {.arg newdata} must be a {.cls data.frame} object."
#'     )
#'     stopifnot_(
#'       !is.null(newdata[[id]]), 
#'       "Can't find grouping variable {.var {id}} in {.arg newdata}."
#'     )
#'     stopifnot_(
#'       !is.null(newdata[[time]]), 
#'       "Can't find time index variable {.var {time}} in {.arg newdata}."
#'     )
#'     object <- update(object, newdata = newdata)
#'   }
#'   S <- object$n_states
#'   M <- object$n_symbols
#'   C <- object$n_channels
#'   N <- object$n_sequences
#'   T_ <- object$length_of_sequences
#'   initial_probs <- get_pi(object$coefficients$gamma_pi_raw, object$X_initial, 0)
#'   transition_probs <- get_A(object$coefficients$gamma_A_raw, object$X_transition, 0)
#'   emission_probs <- if (C == 1) {
#'     get_B(object$coefficients$gamma_B_raw, object$X_emission, 0, 0) 
#'   } else {
#'     get_multichannel_B(object$object$gamma_B_raw, object$X_emission, S, M, 0, 0) 
#'   } 
#'   if (C == 1) {
#'     ids <- rownames(object$observations)
#'     times <- colnames(object$observations)
#'     symbol_names <- list(object$symbol_names)
#'   } else {
#'     ids <- rownames(object$observations[[1]])
#'     times <- colnames(object$observations[[1]])
#'     symbol_names <- object$symbol_names
#'   }
#'   out <- list()
#'   out$initial_probs <- data.frame(
#'     id = rep(ids, each = S),
#'     state = object$state_names,
#'     estimate = c(initial_probs)
#'   )
#'   
#'   out$transition_probs <- data.frame(
#'     id = rep(ids, each = S^2 * T_),
#'     time = rep(times, each = S^2),
#'     state_from = object$state_names,
#'     state_to = rep(object$state_names, each = S),
#'     estimate = unlist(transition_probs)
#'   )
#'   out$emission_probs <- do.call(
#'     rbind, 
#'     lapply(seq_len(C), function(i) {
#'       data.frame(
#'         id = rep(ids, each = S * M[i] * T_),
#'         time = rep(times, each = S * M[i]),
#'         state = object$state_names,
#'         channel = object$channel_names[i],
#'         observation = rep(symbol_names[[i]], each = S),
#'         estimate = unlist(emission_probs[((i - 1) * N + 1):(i * N)])
#'       )
#'     })
#'   )
#'   if (C == 1) out$emission_probs$channel <- NULL
#'   
#'   if (!dontchange_colnames) {
#'     colnames(out$initial_probs)[1] <- object$id_variable
#'     colnames(out$transition_probs)[1] <- object$id_variable
#'     colnames(out$transition_probs)[2] <- object$time_variable
#'     colnames(out$emission_probs)[1] <- object$id_variable
#'     colnames(out$emission_probs)[2] <- object$time_variable
#'   }
#'   out
#' }
#' #' @noRd
#' predict.mnhmm <- function(
#'     object, newdata = NULL, dontchange_colnames = FALSE, ...) {
#'   
#'   if (!is.null(newdata)) {
#'     time <- object$time_variable
#'     id <- object$id_variable
#'     stopifnot_(
#'       is.data.frame(newdata), 
#'       "Argument {.arg newdata} must be a {.cls data.frame} object."
#'     )
#'     stopifnot_(
#'       !is.null(newdata[[id]]), 
#'       "Can't find grouping variable {.var {id}} in {.arg newdata}."
#'     )
#'     stopifnot_(
#'       !is.null(newdata[[time]]), 
#'       "Can't find time index variable {.var {time}} in {.arg newdata}."
#'     )
#'     object <- update(object, newdata = newdata)
#'   }
#'   T_ <- object$length_of_sequences
#'   N <- object$n_sequences
#'   S <- object$n_states
#'   M <- object$n_symbols
#'   C <- object$n_channels
#'   D <- object$n_clusters
#'   gamma_omega_raw <- object$coefficients$gamma_omega_raw
#'   initial_probs <- vector("list", D)
#'   transition_probs <- vector("list", D)
#'   emission_probs <- vector("list", D)
#'   for (d in seq_len(D)) {
#'     gamma_pi_raw <- coef_to_cpp_initial(
#'       matrix(
#'         object$coefficients$gamma_pi_raw[d, ,], 
#'         S - 1, nrow(object$X_initial)
#'       )
#'     )
#'     gamma_A_raw <- coef_to_cpp_transition(
#'       array(
#'         object$coefficients$gamma_A_raw[d, , ,], 
#'         dim = c(S, S - 1, nrow(object$X_transition))
#'       )
#'     )
#'     gamma_B_raw <- coef_to_cpp_emission(
#'       if (C == 1) {
#'         array(
#'           object$coefficients$gamma_B_raw[d, , ,],
#'           dim = c(S, M - 1, nrow(object$X_emission))
#'         )
#'       } else {
#'         object$coefficients$gamma_B_raw[d, ]
#'       },
#'       1,
#'       C > 1
#'     )
#'     initial_probs[[d]] <- get_pi(object$coefficients$gamma_pi_raw, object$X_initial, 0)
#'     transition_probs[[d]] <- get_A(object$coefficients$gamma_A_raw, object$X_transition, 0)
#'     emission_probs[[d]] <- if (C == 1) {
#'       get_B(object$coefficients$gamma_B_raw, object$X_emission, 0, 0) 
#'     } else {
#'       get_multichannel_B(object$coefficients$gamma_B_raw, object$X_emission, S, M, 0, 0) 
#'     }
#'   }
#'   if (C == 1) {
#'     ids <- rownames(object$observations)
#'     times <- colnames(object$observations)
#'     symbol_names <- list(object$symbol_names)
#'   } else {
#'     ids <- rownames(object$observations[[1]])
#'     times <- colnames(object$observations[[1]])
#'     symbol_names <- object$symbol_names
#'   }
#'   out <- list()
#'   out$initial_probs <- data.frame(
#'     cluster = rep(object$cluster_names, each = S * N),
#'     id = rep(ids, each = S),
#'     state = object$state_names,
#'     estimate = unlist(initial_probs)
#'   )
#'   out$transition_probs <- data.frame(
#'     cluster = rep(object$cluster_names, each = S^2 * T_ * N),
#'     id = rep(ids, each = S^2 * T_),
#'     time = rep(times, each = S^2),
#'     state_from = object$state_names,
#'     state_to = rep(object$state_names, each = S),
#'     estimate = unlist(transition_probs)
#'   )
#'   out$emission_probs <- data.frame(
#'     cluster = rep(object$cluster_names, each = S * sum(M) * T_ * N),
#'     id = unlist(lapply(seq_len(C), function(i) rep(ids, each = S * M[i] * T_))),
#'     time = unlist(lapply(seq_len(C), function(i) rep(times, each = S * M[i]))),
#'     state = object$state_names,
#'     channel = rep(object$channel_names, S * M * T_ * N),
#'     observation = rep(unlist(symbol_names), each = S),
#'     estimate = unlist(emission_probs)
#'   )
#'   if (C == 1) emission_probs$channel <- NULL
#'   out$cluster_probs <- data.frame(
#'     cluster = object$cluster_names,
#'     id = rep(ids, each = D),
#'     estimate = c(get_omega(object$coefficients$gamma_omega_raw, 
#'                            object$coefficients$X_cluster, 0))
#'   )
#'   if (!dontchange_colnames) {
#'     colnames(out$initial_probs)[2] <- object$id_variable
#'     colnames(out$transition_probs)[2] <- object$id_variable
#'     colnames(out$transition_probs)[3] <- object$time_variable
#'     colnames(out$emission_probs)[2] <- object$id_variable
#'     colnames(out$emission_probs)[3] <- object$time_variable
#'     colnames(out$cluster_probs)[2] <- object$id_variable
#'   }
#'   out
#' }
