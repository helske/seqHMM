#' Posterior Probabilities forHidden Markov Models
#'
#' Function `posterior_probs` computes the posterior probabilities of hidden 
#' states of a (mixture) hidden Markov model.
#'
#' @export
#' @param modelA hidden Markov model.
#' @param log_space Internally compute posterior probabilities in logarithmic 
#' scale. The default is `TRUE`, which is also only option for 
#' non-homogenous models.
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
posterior_probs.hmm <- function(model, log_space = TRUE, ...) {
  fb <- forward_backward(model, log_space = log_space)
  if (!log_space) {
    fb$forward_probs * fb$backward_probs /
      aperm(replicate(sum(model$n_states), fb$scaling_factors), c(3, 1, 2))
  } else {
    ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
    exp(fb$forward_probs + fb$backward_probs -
          array(
            rep(ll, each = sum(model$n_states) * model$length_of_sequences),
            c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
          ))
  }
}
#' @rdname posterior_probs
#' @export
posterior_probs.mhmm <- function(model, log_space = TRUE, ...) {
  fb <- forward_backward(model, log_space = log_space)
  if (!log_space) {
    fb$forward_probs * fb$backward_probs /
      aperm(replicate(sum(model$n_states), fb$scaling_factors), c(3, 1, 2))
  } else {
    ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
    exp(fb$forward_probs + fb$backward_probs -
          array(
            rep(ll, each = sum(model$n_states) * model$length_of_sequences),
            c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
          )
    )
  }
}
#' @rdname posterior_probs
#' @export
posterior_probs.nhmm <- function(model, ...) {
  fb <- forward_backward(model)
  ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
  exp(fb$forward_probs + fb$backward_probs -
        array(
          rep(ll, each = sum(model$n_states) * model$length_of_sequences),
          c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
        )
  )
}
#' @rdname posterior_probs
#' @export
posterior_probs.mnhmm <- function(model, ...) {
  fb <- forward_backward(model)
  ll <- apply(fb$forward_probs[, model$length_of_sequences,], 2, logSumExp)
  exp(fb$forward_probs + fb$backward_probs -
        array(
          rep(ll, each = sum(model$n_states) * model$length_of_sequences),
          c(sum(model$n_states), model$length_of_sequences, model$n_sequences)
        )
  )
}