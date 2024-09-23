#' Most Probable Paths of Hidden States
#'
#' Function `hidden_paths` computes the most probable path of
#' hidden states of a (mixture) hidden Markov model given the observed sequences.
#'
#' @export
#' @param model A hidden Markov model.
#' @param respect_void If `TRUE` (default), states at the time points
#' corresponding to `TraMineR`'s void in the observed sequences are set to void
#' in the hidden state sequences as well.
#' @param ... Ignored.
#' @return The most probable paths of hidden states as an `stslist` object
#' (see [seqdef()]). The log-probability is included as an attribute 
#' `log_prop`.
#'
#' @examples
#' # Load a pre-defined HMM
#' data("hmm_biofam")
#'
#' # Compute the most probable hidden state paths given the data and the model
#' mpp <- hidden_paths(hmm_biofam)
#'
#' # Plot hidden paths for the first 100 individuals
#' stacked_sequence_plot(mpp, type = "i", ids = 1:100)
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
hidden_paths.hmm <- function(model, respect_void = TRUE, ...) {
  model$initial_probs <- log(model$initial_probs)
  model$transition_probs <- log(model$transition_probs)
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  out <- viterbi(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray
  )
  create_mpp_seq(out, model, respect_void)
}
#' @rdname hidden_paths
#' @export
hidden_paths.mhmm <- function(model, respect_void = TRUE, ...) {
  model <- .combine_models(model)
  model$initial_probs <- log(model$initial_probs)
  model$transition_probs <- log(model$transition_probs)
  obsArray <- create_obsArray(model)
  emissionArray <- create_emissionArray(model)
  out <- viterbix(
    model$transition_probs, emissionArray,
    model$initial_probs, obsArray, model$coefficients,
    model$X, model$n_states_in_clusters
  )
  create_mpp_seq(out, model, respect_void)
}
#' @rdname hidden_paths
#' @export
hidden_paths.nhmm <- function(model, respect_void = TRUE, ...) {
  
  obsArray <- create_obsArray(model)
  if (model$n_channels == 1) {
    out <- viterbi_nhmm_singlechannel(
      model$coefficients$gamma_pi_raw, model$X_initial,
      model$coefficients$gamma_A_raw, model$X_transition,
      model$coefficients$gamma_B_raw, model$X_emission,
      obsArray[1, , ])
  } else {
    out <- viterbi_nhmm_multichannel(
      model$coefficients$gamma_pi_raw, model$X_initial,
      model$coefficients$gamma_A_raw, model$X_transition,
      model$coefficients$gamma_B_raw, model$X_emission,
      obsArray, model$n_symbols)
  }
  create_mpp_seq(out, model, respect_void)
}
#' @rdname hidden_paths
#' @export
hidden_paths.mnhmm <- function(model, respect_void = TRUE, ...) {
  
  obsArray <- create_obsArray(model)
  if (model$n_channels == 1) {
    out <- viterbi_mnhmm_singlechannel(
      model$coefficients$gamma_pi_raw, model$X_initial,
      model$coefficients$gamma_A_raw, model$X_transition,
      model$coefficients$gamma_B_raw, model$X_emission,
      model$coefficients$gamma_omega_raw, model$X_cluster,
      array(obsArray, dim(obsArray)[2:3]))
  } else {
    out <- viterbi_mnhmm_multichannel(
      model$coefficients$gamma_pi_raw, model$X_initial,
      model$coefficients$gamma_A_raw, model$X_transition,
      model$coefficients$gamma_B_raw, model$X_emission,
      model$coefficients$gamma_omega_raw, model$X_cluster,
      obsArray, model$n_symbols)
  }
  if (identical(model$state_names[[1]], model$state_names[[2]])) {
    model$state_names <- paste0(
      rep(model$cluster_names, each = model$n_states), ": ",
      unlist(model$state_names)
    )
  } else {
    model$state_names <- unlist(model$state_names)
  }
  model$n_states <- length(model$state_names)
  create_mpp_seq(out, model, respect_void)
}
#' Create a seqdef Object from the Viterbi Algorithm Output
#' @noRd
create_mpp_seq <- function(out, model, respect_void) {
  if (model$n_sequences == 1) {
    mpp <- model$state_names[out$q + 1]
  } else {
    mpp <- t(apply(out$q + 1, 2, function(x) model$state_names[x]))
  }
  if (model$n_channels == 1) model$observations <- list(model$observations)
  if (respect_void) {
    void_symbol <- attr(model$observations[[1]], "void")
    voids <- vector("list", model$n_channels)
    for (i in 1:model$n_channels) {
      voids[[i]] <- which(model$observations[[i]] == void_symbol)
    }
    mpp[unique(unlist(voids))] <- NA
  }
  mpp <- suppressWarnings(
    suppressMessages(
      seqdef(
        mpp,
        alphabet = model$state_names, id = rownames(model$observations[[1]]),
        start = attr(model$observations[[1]], "start"),
        xtstep = attr(model$observations[[1]], "xtstep"),
        void = void_symbol
      )
    )
  )
  if (sum(model$n_states) <= 200) {
    TraMineR::cpal(mpp) <- seqHMM::colorpalette[[sum(model$n_states)]]
  } else {
    cp <- NULL
    k <- 200
    p <- 0
    while (sum(model$n_states) - p > 0) {
      cp <- c(cp, seqHMM::colorpalette[[k]])
      p <- p + k
      k <- k - 1
    }
    TraMineR::cpal(mpp) <- cp[1:sum(model$n_states)]
  }
  attr(mpp, "log_prob") <- c(out$logp)
  mpp
}
