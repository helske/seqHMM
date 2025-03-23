#' Build a Latent Class Model
#'
#' Function `build_lcm` is a shortcut for constructing a latent class model
#' as a restricted case of an `mhmm` object.
#'
#' @export
#' @param observations An `stslist` object (see [TraMineR::seqdef()]) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param n_clusters A scalar giving the number of clusters/submodels
#' (not used if starting values for model parameters are given with `emission_probs`).
#' @param emission_probs A matrix containing emission probabilities for each class by rows,
#'   or in case of multichannel data a list of such matrices.
#'   Note that the matrices must have dimensions k x s where k is the number of
#'   latent classes and s is the number of unique symbols (observed states) in the
#'   data. Emission probabilities should follow the ordering of the alphabet of
#'   observations (`alphabet(observations)`, returned as `symbol_names`).
#' @param formula Optional formula of class [formula()] for the
#' mixture probabilities. Left side omitted.
#' @param data A data frame containing the variables used in the formula.
#' Ignored if no formula is provided.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the classes/clusters.
#' @param channel_names A vector of optional names for the channels.
#' @return Object of class `mhmm` with the following elements:
#' * `observations`\cr State sequence object or a list of such containing the data.
#' * `transition_probs`\cr A matrix of transition probabilities.
#' * `emission_probs`\cr A matrix or a list of matrices of emission probabilities.
#' * `initial_probs`\cr A vector of initial probabilities.
#' * `coefficients`\cr A matrix of parameter coefficients for covariates (covariates in rows, clusters in columns).
#' * `X`\cr Covariate values for each subject.
#' * `cluster_names`\cr Names for clusters.
#' * `state_names`\cr Names for hidden states.
#' * `symbol_names`\cr Names for observed states.
#' * `channel_names`\cr Names for channels of sequence data
#' * `length_of_sequences`\cr (Maximum) length of sequences.
#' * `sequence_lengths`\cr A vector of sequence lengths.
#' * `n_sequences`\cr Number of sequences.
#' * `n_symbols`\cr Number of observed states (in each channel).
#' * `n_states`\cr Number of hidden states.
#' * `n_channels`\cr Number of channels.
#' * `n_covariates`\cr Number of covariates.
#' * `n_clusters`\cr Number of clusters.
#' 
#' @seealso [fit_model()] for estimating model parameters;
#' [summary.mhmm()] for a summary of a mixture model;
#' [separate_mhmm()] for organizing an `mhmm` object into a list of
#' separate `hmm` objects; and [plot.mhmm()] for plotting
#' mixture models.
#'
#' @examples
#' # Simulate observations from two classes
#' set.seed(123)
#' obs <- seqdef(rbind(
#'   matrix(sample(letters[1:3], 500, TRUE, prob = c(0.1, 0.6, 0.3)), 50, 10),
#'   matrix(sample(letters[1:3], 200, TRUE, prob = c(0.4, 0.4, 0.2)), 20, 10)
#' ))
#'
#' # Initialize the model
#' set.seed(9087)
#' model <- build_lcm(obs, n_clusters = 2)
#'
#' # Estimate model parameters
#' fit <- fit_model(model)
#'
#' # How many of the observations were correctly classified:
#' sum(summary(fit$model)$most_probable_cluster == rep(c("Class 2", "Class 1"), 
#'   times = c(500, 200)))
#'
#' ############################################################
#' \dontrun{
#' # LCM for longitudinal data
#'
#' # Define sequence data
#' data("mvad", package = "TraMineR")
#' mvad_alphabet <- c(
#'   "employment", "FE", "HE", "joblessness", "school",
#'   "training"
#' )
#' mvad_labels <- c(
#'   "employment", "further education", "higher education",
#'   "joblessness", "school", "training"
#' )
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 15:86,
#'   alphabet = mvad_alphabet, states = mvad_scodes,
#'   labels = mvad_labels, xtstep = 6
#' )
#'
#' # Initialize the LCM with random starting values
#' set.seed(7654)
#' init_lcm_mvad1 <- build_lcm(
#'   observations = mvad_seq,
#'   n_clusters = 2, formula = ~male, data = mvad
#' )
#'
#' # Own starting values for emission probabilities
#' emiss <- rbind(rep(1 / 6, 6), rep(1 / 6, 6))
#'
#' # LCM with own starting values
#' init_lcm_mvad2 <- build_lcm(
#'   observations = mvad_seq,
#'   emission_probs = emiss, formula = ~male, data = mvad
#' )
#'
#' # Estimate model parameters (EM algorithm with random restarts)
#' lcm_mvad <- fit_model(init_lcm_mvad1,
#'   control_em = list(restart = list(times = 5))
#' )$model
#'
#' # Plot the LCM
#' plot(lcm_mvad, interactive = FALSE, ncol = 2)
#'
#' ###################################################################
#'
#' # Binomial regression (comparison to glm)
#'
#' require("MASS")
#' data("birthwt")
#'
#' model <- build_lcm(
#'   observations = seqdef(birthwt$low), emission_probs = diag(2),
#'   formula = ~ age + lwt + smoke + ht, data = birthwt
#' )
#' fit <- fit_model(model)
#' summary(fit$model)
#' summary(glm(low ~ age + lwt + smoke + ht, binomial, data = birthwt))
#'
#'
#' # Multinomial regression (comparison to multinom)
#'
#' require("nnet")
#'
#' set.seed(123)
#' n <- 100
#' X <- cbind(1, x1 = runif(n, 0, 1), x2 = runif(n, 0, 1))
#' coefs <- cbind(0, c(-2, 5, -2), c(0, -2, 2))
#' pr <- exp(X %*% coefs) + stats::rnorm(n * 3)
#' pr <- pr / rowSums(pr)
#' y <- apply(pr, 1, which.max)
#' table(y)
#'
#' model <- build_lcm(
#'   observations = seqdef(y), emission_probs = diag(3),
#'   formula = ~ x1 + x2, data = data.frame(X[, -1])
#' )
#' fit <- fit_model(model)
#' summary(fit$model)
#' summary(multinom(y ~ x1 + x2, data = data.frame(X[, -1])))
#' }
build_lcm <- function(observations, n_clusters, emission_probs,
                      formula = NULL, data = NULL, coefficients = NULL,
                      cluster_names = NULL, channel_names = NULL) {
  
  n_clusters_given <- !missing(n_clusters)
  emission_probs_given <- !missing(emission_probs)
  stopifnot_(
    n_clusters_given || emission_probs_given,
    "Provide either {.arg emission_probs} or {.arg n_clusters}."
  )
  observations <- .check_observations(observations, channel_names)
  n_channels <- attr(observations, "n_channels")
  n_symbols <- attr(observations, "n_symbols")
  channel_names <- attr(observations, "channel_names")
  
  if (!emission_probs_given) {
    emission_probs <- simulate_emission_probs(1L, n_symbols, n_clusters)
  } else {
    if (n_channels > 1L) {
      n_clusters <- nrow(emission_probs[[1]])
      emission_probs_list <- vector("list", n_clusters)
      for (i in seq_len(n_channels)) {
        stopifnot_(
          nrow(emission_probs[[i]]) != n_clusters,
          "Different number of rows in the list components of {.arg emission_probs}."
        )
      }
      for (i in seq_len(n_clusters)) {
        emission_probs_list[[i]] <- vector("list", n_channels)
        for (j in seq_len(n_channels)) {
          emission_probs_list[[i]][[j]] <- emission_probs[[j]][i, , drop = FALSE]
        }
      }
    } else {
      n_clusters <- nrow(emission_probs)
      emission_probs_list <- vector("list", n_clusters)
      for (i in 1:n_clusters) {
        emission_probs_list[[i]] <- emission_probs[i, , drop = FALSE]
      }
    }
    emission_probs <- emission_probs_list
  }
  if (is.null(cluster_names)) {
    cluster_names <- paste("Class", seq_len(n_clusters))
  } else if (length(cluster_names) != n_clusters) {
    warning_(
      "The length of {.arg cluster_names} does not match the number of clusters.
      Names were not used."
      )
    cluster_names <- paste("Class", seq_len(n_clusters))
  }
  n_states <- rep(1, n_clusters)
  initial_probs <- replicate(n_clusters, 1, simplify = FALSE)
  transition_probs <- replicate(n_clusters, matrix(1), simplify = FALSE)

  model <- build_mhmm(
    observations = observations, transition_probs = transition_probs,
    emission_probs = emission_probs, initial_probs = initial_probs,
    formula = formula, data = data, coefficients = coefficients,
    cluster_names = cluster_names, state_names = cluster_names,
    channel_names = channel_names)
  model$call <- match.call()
  attr(model, "type") <- "lcm"
  model
}
