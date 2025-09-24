#' Print Method for a Hidden Markov Model
#'
#' Prints the parameters of a (mixture) hidden Markov model.
#'
#' @export
#' @rdname print
#' @param x Hidden Markov model.
#' @param digits Minimum number of significant digits to print.
#' @param ... Further arguments to `print.default`.
#' @seealso [build_hmm()] and [fit_model()] for building and
#'   fitting hidden Markov models.
print.hmm <- function(x, digits = 3, ...) {
  if (x$n_channels == 1) {
    if (attr(x, "type") == "mm") {
      cat("\nMarkov model: \n")
      cat("\nNumber of sequences:", x$n_sequences)
      cat("\nNumber of time points:", x$length_of_sequences)
      cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
      print.listof(list(
        "\nInitial probabilities" = x$initial_probs,
        "Transition probabilities" = x$transition_probs
      ), digits = digits, ...)
    } else {
      cat("\nHidden Markov model: \n")
      cat("\nNumber of sequences:", x$n_sequences)
      cat("\nNumber of time points:", x$length_of_sequences)
      cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
      cat("\nNumber of hidden states:", x$n_states, "\n")
      print.listof(list(
        "Initial probabilities" = x$initial_probs,
        "Transition probabilities" = x$transition_probs,
        "Emission probabilities" = x$emission_probs
      ), digits = digits, ...)
    }
  } else {
    cat("\nHidden Markov model: \n")
    cat("\nNumber of sequences:", x$n_sequences)
    cat("\nNumber of time points:", x$length_of_sequences)
    cat("\nNumber of observation channels:", x$n_channels)
    cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
    cat("\nNumber of hidden states:", x$n_states, "\n")
    cat("\n")
    print.listof(list("Initial probabilities" = x$initial_probs), digits = digits, ...)
    cat("\n")
    print.listof(list("Transition probabilities" = x$transition_probs), digits = digits, ...)
    cat("\n")
    cat("Emission probabilities :\n")
    print.listof(x$emission_probs, digits = digits, ...)
    cat("\n")
  }
}
#' @export
#' @rdname print
print.mhmm <- function(x, digits = 3, ...) {
  if (attr(x, "type") == "lcm") {
    cat("Latent class model: \n")
  } else {
    if (attr(x, "type") == "mmm") {
      cat("Mixture Markov model: \n")
    } else {
      cat("Mixture hidden Markov model: \n")
    }
  }
  cat("\nNumber of sequences:", x$n_sequences)
  cat("\nNumber of time points:", x$length_of_sequences)
  if (attr(x, "type") == "lcm") {
    cat("\nNumber of observation channels:", x$n_channels)
    cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
    cat("\nNumber of clusters:", x$n_clusters)
  } else {
    if (attr(x, "type") == "mmm") {
      cat("\nNumber of states:", x$n_states)
      cat("\nNumber of clusters:", x$n_clusters)
    } else {
      cat("\nNumber of observation channels:", x$n_channels)
      cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
      cat("\nNumber of states:", paste(x$n_states, collapse = ", "))
      cat("\nNumber of clusters:", x$n_clusters)
    }
  }
  cat("\nCoefficients :\n")
  print(x$coefficients, digits = digits, ...)
  
  if (attr(x, "type") != "lcm") {
    cat("\nInitial probabilities :\n")
    print.listof(x$initial_probs, digits = digits, ...)
    cat("Transition probabilities :\n")
    print.listof(x$transition_probs, digits = digits, ...)
  } else {
    cat("\n")
  }
  if (attr(x, "type") != "mmm") {
    cat("Emission probabilities :\n")
    
    if (x$n_channels == 1) {
      if (attr(x, "type") == "lcm") {
        tmp <- do.call(rbind, x$emission_probs)
        rownames(tmp) <- x$state_names
        print(tmp, digits = digits, ...)
      } else {
        print.listof(x$emission_probs, digits = digits, ...)
      }
    } else {
      if (attr(x, "type") == "lcm") {
        for (i in 1:x$n_channels) {
          cat(x$channel_names[i], ":\n")
          print(do.call(rbind, sapply(x$emission_probs, "[", i)), digits = digits, ...)
          cat("\n")
        }
      } else {
        for (i in 1:x$n_clusters) {
          cat(names(x$emission_probs)[i], ":\n\n")
          print.listof(x$emission_probs[[i]], digits = digits, ...)
        }
      }
    }
  }
  
  cat("\n")
}
#' @export
#' @rdname print
print.nhmm <- function(x, digits = 3, ...) {
  fanhmm <- inherits(x, "fanhmm")
  if (fanhmm) {
    cat("Non-homogeneous feedback-augmented hidden Markov model: \n")
  } else {
    cat("Non-homogeneous hidden Markov model: \n")
  }
  cat("\nNumber of sequences:", x$n_sequences)
  if (fanhmm && identical(x$prior_obs, 0L)) {
    cat("\nNumber of time points:", x$length_of_sequences, " (after fixing the first time point)")
  } else {
    cat("\nNumber of time points:", x$length_of_sequences)
  }
  cat("\nNumber of observation channels:", x$n_channels)
  cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
  cat("\nNumber of hidden states:", x$n_states)
  
  cat("\nFormula for initial state probabilities: \n")
  print(x$initial_formula, showEnv = FALSE)
  cat("\nFormula for transition probabilities: \n")
  print(x$transition_formula, showEnv = FALSE)
  cat("\nFormula(s) for emission probabilities: \n")
  for (y in x$responses) {
    print(x$emission_formula[[y]], showEnv = FALSE)
  }
}
#' @export
#' @rdname print
print.mnhmm <- function(x, digits = 3, ...) {
  fanhmm <- inherits(x, "fanhmm")
  if (fanhmm) {
    cat("Mixture of non-homogeneous feedback-augmented hidden Markov models: \n")
  } else {
    cat("Mixture of non-homogeneous hidden Markov models: \n")
  }
  cat("\nNumber of sequences:", x$n_sequences)
  if (fanhmm && identical(x$prior_obs, 0L)) {
    cat("\nNumber of time points:", x$length_of_sequences, " (after fixing the first time point)")
  } else {
    cat("\nNumber of time points:", x$length_of_sequences)
  }
  cat("\nNumber of observation channels:", x$n_channels)
  cat("\nNumber of observed symbols:", paste(x$n_symbols, collapse = ", "))
  cat("\nNumber of hidden states:", x$n_states)
  cat("\nNumber of clusters:", x$n_clusters)
  
  cat("\nFormula for initial state probabilities: \n")
  print(x$initial_formula, showEnv = FALSE)
  cat("\nFormula for transition probabilities: \n")
  print(x$transition_formula, showEnv = FALSE)
  cat("\nFormula(s) for emission probabilities: \n")
  for (y in x$responses) {
    print(x$emission_formula[[y]], showEnv = FALSE)
  }
  cat("\nFormula for mixture probabilities: \n")
  print(x$cluster_formula, showEnv = FALSE)
}
