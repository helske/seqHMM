#' Posterior Probabilities forHidden Markov Models
#'
#' Function `posterior_probs` computes the posterior probabilities of hidden 
#' states of a (mixture) hidden Markov model.
#'
#' @export
#' @param model A hidden Markov model.
#' @param log_space Internally compute posterior probabilities in logarithmic 
#' scale. The default is `TRUE`, which is also only option for 
#' non-homogenous models.
#' @param as_data_frame If `TRUE` (default), the output is returned as a 
#' data.frame. Otherwise, a 3d array is returned.
#' @param ... Ignored.
#' @return Posterior probabilities. In case of multiple observations,
#' these are computed independently for each sequence.
#' @examples
#' # Load a pre-defined MHMM
#' data("mhmm_biofam")
#'
#' # Compute posterior probabilities
#' pb <- posterior_probs(mhmm_biofam)
#'
#' # Locally most probable states for the first subject:
#' pb[, , 1]
#' @rdname posterior_probs
#' @export
posterior_probs <- function(model, ...) {
  UseMethod("posterior_probs", model)
}
#' @rdname posterior_probs
#' @export
posterior_probs.hmm <- function(model, log_space = TRUE, as_data_frame = TRUE, 
                                ...) {
  fb <- forward_backward(model, log_space = log_space, as_data_frame = FALSE)
  if (!log_space) {
    out <- fb$forward_probs * fb$backward_probs /
      aperm(replicate(sum(model$n_states), fb$scaling_factors), c(3, 1, 2))
  } else {
    ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
    out <- exp(fb$forward_probs + fb$backward_probs -
                 array(
                   rep(ll, each = sum(model$n_states) * model$length_of_sequences),
                   c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
                 ))
  }
  if (as_data_frame) {
    cbind(expand.grid(dimnames(out)), probability = c(out))  
  } else {
    out
  }
}
#' @rdname posterior_probs
#' @export
posterior_probs.mhmm <- function(model, log_space = TRUE, as_data_frame = TRUE, 
                                 ...) {
  fb <- forward_backward(model, log_space = log_space, as_data_frame = FALSE)
  if (!log_space) {
    out <- fb$forward_probs * fb$backward_probs /
      aperm(replicate(sum(model$n_states), fb$scaling_factors), c(3, 1, 2))
  } else {
    ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
    out <- exp(fb$forward_probs + fb$backward_probs -
                 array(
                   rep(ll, each = sum(model$n_states) * model$length_of_sequences),
                   c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
                 )
    )
  }
  if (as_data_frame) {
    cbind(expand.grid(dimnames(out)), probability = c(out))  
  } else {
    out
  }
}
#' @rdname posterior_probs
#' @export
posterior_probs.nhmm <- function(model, as_data_frame = TRUE, ...) {
  fb <- forward_backward(model, as_data_frame = FALSE)
  ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
  out <- exp(fb$forward_probs + fb$backward_probs -
               array(
                 rep(ll, each = model$n_states * model$length_of_sequences),
                 c(model$n_states, model$length_of_sequences, model$n_sequences)
               )
  )
  if (as_data_frame) {
    cbind(expand.grid(dimnames(out)), probability = c(out))  
  } else {
    out
  }
}
#' @rdname posterior_probs
#' @export
posterior_probs.mnhmm <- function(model, as_data_frame = TRUE,...) {
  D <-  model$n_clusters
  S <- model$n_states
  T <- model$length_of_sequences
  fb <- forward_backward(model, as_data_frame = FALSE)
  ll <- apply(fb$forward_probs[, T, ], 2, logSumExp)
  out <- exp(
    fb$forward_probs + fb$backward_probs -
      array(rep(ll, each = D * S * T), c(D * S, T, model$n_sequences))
  )
  if (as_data_frame) {
    cbind(expand.grid(dimnames(out)), probability = c(out))  
  } else {
    out
  }
}