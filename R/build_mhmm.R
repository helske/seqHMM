#' Build a Mixture Hidden Markov Model
#'
#' Function \code{build_mhmm} constructs a mixture hidden Markov model object of class \code{mhmm}.
#'
#' The returned model contains some attributes such as \code{nobs} and \code{df},
#' which define the number of observations in the  model and the number of estimable
#' model parameters, used in computing BIC.
#' When computing \code{nobs} for a multichannel model with \eqn{C} channels,
#' each observed value in a single channel amounts to \eqn{1/C} observation,
#' i.e. a fully observed time point for a single sequence amounts to one observation.
#' For the degrees of freedom \code{df}, zero probabilities of the initial model are
#' defined as structural zeroes.
#'
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param n_states A numerical vector giving the number of hidden states in each submodel
#' (not used if starting values for model parameters are given with
#' \code{initial_probs}, \code{transition_probs}, or \code{emission_probs}).
#' @param transition_probs A list of matrices of transition
#'   probabilities for the submodel of each cluster.
#' @param emission_probs A list which contains matrices of emission probabilities or
#'   a list of such objects (one for each channel) for the submodel of each cluster.
#'   Note that the matrices must have dimensions \eqn{m x s} where \eqn{m} is the number of
#'   hidden states and \eqn{s} is the number of unique symbols (observed states) in the
#'   data. Emission probabilities should follow the ordering of the alphabet of
#'   observations (\code{alphabet(observations)}, returned as \code{symbol_names}).
#' @param initial_probs A list which contains vectors of initial state
#'   probabilities for the submodel of each cluster.
#' @param formula Optional formula of class \code{\link{formula}} for the
#' mixture probabilities. Left side omitted.
#' @param data A data frame containing the variables used in the formula.
#' Ignored if no formula is provided.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the clusters.
#' @param state_names A list of optional labels for the hidden states. If \code{NULL},
#' the state names are taken as row names of transition matrices. If this is also \code{NULL},
#' numbered states are used.
#' @param channel_names A vector of optional names for the channels.
#' @param ... Additional arguments to \code{simulate_transition_probs}.
#' @return Object of class \code{mhmm} with following elements:
#' \describe{
#'    \item{\code{observations}}{State sequence object or a list of such containing the data.}
#'    \item{\code{transition_probs}}{A matrix of transition probabilities.}
#'    \item{\code{emission_probs}}{A matrix or a list of matrices of emission probabilities.}
#'    \item{\code{initial_probs}}{A vector of initial probabilities.}
#'    \item{\code{coefficients}}{A matrix of parameter coefficients for covariates (covariates in rows, clusters in columns).}
#'    \item{\code{X}}{Covariate values for each subject.}
#'    \item{\code{cluster_names}}{Names for clusters.}
#'    \item{\code{state_names}}{Names for hidden states.}
#'    \item{\code{symbol_names}}{Names for observed states.}
#'    \item{\code{channel_names}}{Names for channels of sequence data}
#'    \item{\code{length_of_sequences}}{(Maximum) length of sequences.}
#'    \item{\code{n_sequences}}{Number of sequences.}
#'    \item{\code{n_symbols}}{Number of observed states (in each channel).}
#'    \item{\code{n_states}}{Number of hidden states.}
#'    \item{\code{n_channels}}{Number of channels.}
#'    \item{\code{n_covariates}}{Number of covariates.}
#'    \item{\code{n_clusters}}{Number of clusters.}
#' }
#' @seealso \code{\link{fit_model}} for fitting mixture Hidden Markov models;
#' \code{\link{summary.mhmm}} for a summary of a MHMM; \code{\link{separate_mhmm}} for
#' reorganizing a MHMM into a list of separate hidden Markov models; and
#' \code{\link{plot.mhmm}} for plotting \code{mhmm} objects.
#'
#' @references Helske S. and Helske J. (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R,
#' Journal of Statistical Software, 88(3), 1-32. doi:10.18637/jss.v088.i03
#'
#' @examples
#'
#' data("biofam3c")
#'
#' ## Building sequence objects
#' marr_seq <- seqdef(biofam3c$married,
#'   start = 15,
#'   alphabet = c("single", "married", "divorced")
#' )
#' child_seq <- seqdef(biofam3c$children,
#'   start = 15,
#'   alphabet = c("childless", "children")
#' )
#' left_seq <- seqdef(biofam3c$left,
#'   start = 15,
#'   alphabet = c("with parents", "left home")
#' )
#'
#' ## Choosing colors
#' attr(marr_seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(child_seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(left_seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#'
#' ## MHMM with random starting values, no covariates
#' set.seed(468)
#' init_mhmm_bf1 <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   n_states = c(4, 4, 6),
#'   channel_names = c("Marriage", "Parenthood", "Residence")
#' )
#'
#'
#' ## Starting values for emission probabilities
#'
#' # Cluster 1
#' B1_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.3, 0.6, 0.1, # High probability for married
#'     0.3, 0.3, 0.4
#'   ), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE
#' )
#'
#' B1_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.9, 0.1
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' B1_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9, # High probability for having left home
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' # Cluster 2
#'
#' B2_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.7, 0.2, 0.1
#'   ),
#'   nrow = 4, ncol = 3, byrow = TRUE
#' )
#'
#' B2_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' B2_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 4, ncol = 2, byrow = TRUE
#' )
#'
#' # Cluster 3
#' B3_marr <- matrix(
#'   c(
#'     0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.3, 0.4, 0.3,
#'     0.1, 0.1, 0.8
#'   ), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE
#' )
#'
#' B3_child <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9
#'   ),
#'   nrow = 6, ncol = 2, byrow = TRUE
#' )
#'
#'
#' B3_left <- matrix(
#'   c(
#'     0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9,
#'     0.1, 0.9
#'   ),
#'   nrow = 6, ncol = 2, byrow = TRUE
#' )
#'
#' # Starting values for transition matrices
#' A1 <- matrix(
#'   c(
#'     0.80, 0.16, 0.03, 0.01,
#'     0, 0.90, 0.07, 0.03,
#'     0, 0, 0.90, 0.10,
#'     0, 0, 0, 1
#'   ),
#'   nrow = 4, ncol = 4, byrow = TRUE
#' )
#'
#' A2 <- matrix(
#'   c(
#'     0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
#'     0, 0.70, 0.10, 0.10, 0.05, 0.05,
#'     0, 0, 0.85, 0.01, 0.10, 0.04,
#'     0, 0, 0, 0.90, 0.05, 0.05,
#'     0, 0, 0, 0, 0.90, 0.10,
#'     0, 0, 0, 0, 0, 1
#'   ),
#'   nrow = 6, ncol = 6, byrow = TRUE
#' )
#'
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#'
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort,
#'   labels = c("1909-1935", "1936-1945", "1946-1957")
#' )
#'
#' ## MHMM with own starting values and covariates
#' init_mhmm_bf2 <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   transition_probs = list(A1, A1, A2),
#'   emission_probs = list(
#'     list(B1_marr, B1_child, B1_left),
#'     list(B2_marr, B2_child, B2_left),
#'     list(B3_marr, B3_child, B3_left)
#'   ),
#'   formula = ~ sex + cohort, data = biofam3c$covariates,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Marriage", "Parenthood", "Residence"),
#'   state_names = list(
#'     paste("State", 1:4), paste("State", 1:4),
#'     paste("State", 1:6)
#'   )
#' )
#'
build_mhmm <- function(observations,
                       n_states,
                       transition_probs,
                       emission_probs,
                       initial_probs,
                       formula = NULL,
                       data = NULL,
                       coefficients = NULL,
                       cluster_names = NULL,
                       state_names = NULL,
                       channel_names = NULL, ...) {
  
  observations <- check_observations(observations, channel_names)
  n_channels <- attr(observations, "n_channels")
  n_symbols <- attr(observations, "n_symbols")
  channel_names <- attr(observations, "channel_names")
  symbol_names <- attr(observations, "symbol_names")
  
  inits_given <- !missing(transition_probs) || !missing(initial_probs) || !missing(emission_probs)
  n_states_given <- !missing(n_states)
  stopifnot_(
    inits_given || n_states_given,
    "Provide either {.arg n_states} or all three of {.arg initial_probs}, ",
    "{.arg transition_probs}, and {.arg emission_probs}."
  )
  
  # if any initial values are given, ignore n_states and use these
  if (inits_given) {
    stopifnot_(
      is.list(transition_probs),
      "{.arg transition_probs} is not a {.cls list}."
    )
    stopifnot_(
      is.list(emission_probs),
      "{.arg emission_probs} is not a {.cls list}."
    )
    stopifnot_(
      is.list(initial_probs),
      "{.arg initial_probs} is not a {.cls list}."
    )
    n_clusters <- length(transition_probs)
    stopifnot_(
      n_clusters > 1L,
      "{.arg transition_probs} is a list of length 1 leading to an
        ordinary HMM. Please use {.fn build_hmm} instead."
    )
    stopifnot_(
      length(emission_probs) == n_clusters,
      "Length of a {.arg emission_probs} does not match with the length of
      {.arg transition_probs}."
    )
    stopifnot_(
      length(initial_probs) == n_clusters,
      "Length of a {.arg initial_probs} does not match with the length of
      {.arg transition_probs}."
    )
    if (is.null(cluster_names)) {
      cluster_names <- paste("Cluster", seq_len(n_clusters))
    } else {
      if (!identical(length(cluster_names), n_clusters)) {
        warning_(
          "The length of {.arg cluster_names} does not match the number of
          clusters. Names were not used."
        )
      }
      cluster_names <- paste("Cluster", seq_len(n_clusters))
    }
    n_states <- integer(n_clusters)
    for (i in seq_len(n_clusters)) {
      transition_probs[[i]] <- 
        check_transition_probs(transition_probs[[i]], state_names[[i]])
      n_states[i] <- nrow(transition_probs[[i]])
      state_names[[i]] <- rownames(transition_probs[[i]])
      initial_probs[[i]] <- 
        check_initial_probs(initial_probs[[i]], n_states[i], state_names[[i]])
      emission_probs <- check_emission_probs(
        emission_probs[[i]], n_states[i], n_channels, n_symbols[i], 
        state_names[[i]], symbol_names[[i]]
      )
    }
  } else {
    # Simulate starting values
    n_clusters <- length(n_states)
    if (identical(n_clusters, 1L)) {
      stop(paste0("Argument 'n_states' is of length 1, leading to ordinary ", 
                  "HMM. Please use 'build_hmm' instead."))
    }
    if (is.null(cluster_names)) {
      cluster_names <- paste("Cluster", seq_len(n_clusters))
    } else if (length(cluster_names) != n_clusters) {
      warning(paste0("The length of argument cluster_names does not match ", 
                     "the number of clusters. Names were not used."))
      cluster_names <- paste("Cluster", seq_len(n_clusters))
    }
    
    transition_probs <- 
      simulate_transition_probs(n_states = n_states, n_clusters = n_clusters,
                                ...)
    initial_probs <- simulate_initial_probs(
      n_states = n_states, 
      n_clusters = n_clusters)
    emission_probs <- simulate_emission_probs(
      n_states = n_states, n_symbols = n_symbols, n_clusters = n_clusters
    )
    for (i in seq_len(n_clusters)) {
      transition_probs[[i]] <- check_transition_probs(transition_probs[[i]], state_names[[i]])
      state_names[[i]] <- rownames(transition_probs[[i]])
      initial_probs[[i]] <- check_initial_probs(initial_probs[[i]], n_states[i], state_names[[i]])
      emission_probs <- check_emission_probs(
        emission_probs[[i]], n_states[i], n_channels, n_symbols[i], state_names[[i]], symbol_names[[i]]
      )
    }
  }
  
  if (is.null(formula)) {
    formula <- stats::formula(~ 1)
    X <- model.matrix(formula, data = data.frame(y = rep(1, n_sequences)))
    n_covariates <- 1L
  } else {
    if (inherits(formula, "formula")) {
      if (is.null(data)) {
        stop("Argument 'data' is missing, but 'formula' was provided.")
      }
      X <- model.matrix(formula, data)
      if (nrow(X) != n_sequences) {
        if (length(all.vars(formula)) > 0 && sum(!complete.cases(data[all.vars(formula)])) > 0) {
          stop(
            paste0(
              "Missing cases are not allowed in covariates. Use e.g. the ",
              "'complete.cases' function to detect them, then fix, impute, or remove."
            )
          )
        } else {
          stop(
            paste0(
              "Number of subjects in data for covariates does not match the ",
              "number of subjects in the sequence data."
            )
          )
        }
      }
      n_covariates <- ncol(X)
    } else {
      stop("Object given for argument formula is not of class formula.")
    }
  }
  if (is.null(coefficients)) {
    coefficients <- matrix(0, n_covariates, n_clusters)
  } else {
    if (ncol(coefficients) != n_clusters | nrow(coefficients) != n_covariates) {
      stop("Wrong dimensions of coefficients.")
    }
    coefficients[, 1] <- 0
  }
  rownames(coefficients) <- colnames(X)
  colnames(coefficients) <- cluster_names
  names(transition_probs) <- names(emission_probs) <- names(initial_probs) <- cluster_names
  
  model <- structure(
    list(
      observations = observations, transition_probs = transition_probs,
      emission_probs = emission_probs, initial_probs = initial_probs,
      coefficients = coefficients, X = X, cluster_names = cluster_names, state_names = state_names,
      symbol_names = symbol_names, channel_names = channel_names,
      length_of_sequences = length_of_sequences,
      n_sequences = n_sequences, n_clusters = n_clusters,
      n_symbols = n_symbols, n_states = n_states,
      n_channels = n_channels,
      n_covariates = n_covariates, formula = formula
    ),
    class = "mhmm",
    nobs = attr(observations, "nobs"),
    df = sum(unlist(initial_probs) > 0) - n_clusters + sum(unlist(transition_probs) > 0) - sum(n_states) +
      sum(unlist(emission_probs) > 0) - sum(n_states) * n_channels + n_covariates * (n_clusters - 1),
    type = "mhmm"
  )
  model
}
