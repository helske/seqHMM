#' Most Probable Paths of Hidden States
#'
#' Function `hidden_paths` computes the most probable path of
#' hidden states of a (mixture) hidden Markov model given the observed sequences.
#'
#' @export
#' @param model A hidden Markov model.
#' @param as_stslist Logical. If `TRUE`, the output the is converted to an 
#' `stslist` object. Default is `FALSE`, which returns a `data.table`.
#' @param ... Ignored.
#' @return The most probable paths of hidden states as an `data.table`. 
#' The log-probability is included as an attribute 
#' `log_prop`.
#'
#' @examples
#' # Load a pre-defined HMM
#' data("hmm_biofam")
#'
#' # Compute the most probable hidden state paths given the data and the model
#' mpp <- hidden_paths(hmm_biofam)
#' head(mpp)
#' # Plot hidden paths for the first 100 individuals
#' seqs <- data_to_stslist(mpp, "id", "time", "state")
#' stacked_sequence_plot(seqs, type = "i", ids = 1:100)
#'
#' # Because the model structure is so sparse that the posterior probabilities are
#' # mostly peaked to single state at each time point, the joint probability of
#' # observations and most probable paths of hidden states is almost identical to
#' # log-likelihood:
#'
#' sum(attr(mpp, "log_prob"))
#' logLik(hmm_biofam)
#'
#' @seealso [hmm_biofam] for information on the model used in the example;
#' and [ggseqplot::ggseqiplot()] and [stacked_sequence_plot()]
#' for plotting hidden paths.
#' 
#' @rdname hidden_paths
#' @export
hidden_paths <- function(model, ...) {
  UseMethod("hidden_paths", model)
}
#' @rdname hidden_paths
#' @export
hidden_paths.hmm <- function(model, as_stslist = FALSE, ...) {
  model$initial_probs <- log(model$initial_probs)
  model$transition_probs <- log(model$transition_probs)
  obsArray <- create_obsArray(model)
  emissionArray <- log(create_emissionArray(model))
  out <- viterbi(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray
  )
  create_mpp_data(out, model, as_stslist)
}
#' @rdname hidden_paths
#' @export
hidden_paths.mhmm <- function(model, as_stslist = FALSE, ...) {
  model <- .combine_models(model)
  model$initial_probs <- log(model$initial_probs)
  model$transition_probs <- log(model$transition_probs)
  obsArray <- create_obsArray(model)
  emissionArray <- log(create_emissionArray(model))
  out <- viterbix(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray, model$coefficients,
    model$X, model$n_states_in_clusters
  )
  create_mpp_data(out, model, as_stslist)
}
#' @rdname hidden_paths
#' @export
hidden_paths.nhmm <- function(model, as_stslist = FALSE, ...) {
  
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    out <- Rcpp_viterbi_fanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B,
      model$prior_obs, model$W_X_B
    )
  } else {
    out <- Rcpp_viterbi_nhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, 
      io(model$X_pi), io(model$X_A), io(model$X_B),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B
    )
  }
  create_mpp_data(out, model, as_stslist)
}
#' @rdname hidden_paths
#' @export
hidden_paths.mnhmm <- function(model, as_stslist = FALSE, ...) {
  
  obs <- create_obs(model)
  if (inherits(model, "fanhmm")) {
    out <- Rcpp_viterbi_mfanhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega, model$prior_obs, model$W_X_B
    )
  } else {
    out <- Rcpp_viterbi_mnhmm(
      obs, model$sequence_lengths, model$n_symbols, 
      model$X_pi, model$X_A, model$X_B, model$X_omega,
      io(model$X_pi), io(model$X_A), io(model$X_B), io(model$X_omega),
      iv(model$X_A), iv(model$X_B), tv(model$X_A), tv(model$X_B),
      model$gammas$gamma_pi, model$gammas$gamma_A, model$gammas$gamma_B, 
      model$gammas$gamma_omega
    )
  }
  model$original_state_names <- model$state_names
  model$state_names <- paste0(
    rep(model$cluster_names, each = model$n_states), ": ",
    unlist(model$state_names)
  )
  create_mpp_data(out, model, as_stslist)
}
#' Create a data.table from the Viterbi Algorithm Output
#' @noRd
create_mpp_data <- function(out, model, as_stslist = FALSE) {
  # avoid CRAN check warnings due to NSE
  state <- cluster <- NULL
  
  if (inherits(model, "nhmm") || inherits(model, "mnhmm")) {
    id <- model$id_variable
    time <- model$time_variable
    d <- model$data[, list(id, time), env = list(id = id, time = time)]
    set(d, j = "state", value = model$state_names[unlist(out$q) + 1L])
  } else {
    id <- "id"
    time <- "time"
    if (model$n_sequences == 1) {
      mpp <- model$state_names[out$q + 1]
    } else {
      mpp <- apply(out$q + 1, 2, \(x) model$state_names[x])
    }
    if (model$n_channels == 1) model$observations <- list(model$observations)
    if (is.null(times <- as.numeric(colnames(model$observations[[1]])))) {
      times <- seq_len(model$length_of_sequences)
    }
    if (is.null(ids <- rownames(model$observations[[1]]))) {
      ids <- seq_len(model$n_sequences)
    }
    ids <- as_factor(ids)
    names(model$sequence_lengths) <- ids
    d <- data.table(
      expand.grid(
        time = times,
        id = ids,
        stringsAsFactors = FALSE
      )[, 2:1],
      state = c(mpp)
    )[time <= times[model$sequence_lengths[id]], ]
    setkey(d, id, time)
  }
  if (as_stslist) {
    d <- suppressMessages(
      data_to_stslist(d, id, time, "state")
    )
  } else {
    if (inherits(model, "combined_mhmm") || inherits(model, "mnhmm")) {
      d[, c("cluster", "state") := tstrsplit(state, " *: *", fixed = FALSE)]
      setcolorder(
        d, 
        c(setdiff(names(d), c("state", "cluster")), c("state", "cluster"))
      )
      d[, cluster := factor(cluster, levels = model$cluster_names)]
      d[, state := factor(state, 
                          levels = unique(unlist(model$original_state_names)))]
    }
  }
  attr(d, "log_prob") <- c(out$logp)
  d
}
