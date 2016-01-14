#' Build a Mixture Markov Model
#'
#' Function \code{build_mmm} is a shortcut for constructing a mixture Markov
#' model as a restricted case of an \code{mhmm} object.
#'
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences.
#' @param transition_probs A list of matrices of transition
#'   probabilities for submodels of each cluster.
#' @param initial_probs A list which contains vectors of initial state
#'   probabilities for submodels of each cluster.
#' @param formula Covariates as an object of class \code{\link{formula}},
#' left side omitted.
#' @param data An optional data frame, list or environment containing the variables
#' in the model. If not found in data, the variables are taken from
#' \code{environment(formula)}.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for 
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number 
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the clusters.
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
#'}
#' @seealso \code{\link{fit_model}} for estimating model parameters;
#' \code{\link{summary.mhmm}} for a summary of a mixture model; 
#' \code{\link{separate_mhmm}} for organizing an \code{mhmm} object into a list of 
#' separate \code{hmm} objects; and \code{\link{plot.mhmm}} for plotting 
#' mixture models.
#'
#' @examples
#'
#'
#' # Define sequence data
#' data("mvad", package = "TraMineR")
#' mvad_alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'   "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'   "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet, states = mvad_scodes,
#'   labels = mvad_labels, xtstep = 6)
#'
#' # Initialize the MMM
#' set.seed(123)
#' mmm_mvad <- build_mmm(observations = mvad_seq,
#'   transition_probs = simulate_transition_probs(n_states = 6, n_clusters = 2),
#'   initial_probs = replicate(2, rep(1/6, 6), simplify = FALSE),
#'   formula = ~male, data = mvad)
#'
#' \dontrun{
#' # Estimate model parameters
#' mmm_mvad <- fit_model(mmm_mvad)$model
#'
#' # Plot model (both clusters in the same plot)
#' require(igraph)
#' plot(mmm_mvad, interactive = FALSE,
#'   # Modify legend position and properties
#'   withlegend = "right", legend.prop = 0.3, cex.legend = 1.2,
#'   # Define vertex layout
#'   layout = layout_as_star,
#'   # Modify edge properties
#'   edge.label = NA, edge.arrow.size = 0.8, edge.curved = 0.2,
#'   # Modify vertex label positions (initial probabilities)
#'   vertex.label.pos = c("left", "right", "right", "left", "left", "right"))
#'   
#' # Summary of the MMM
#' summary(mmm_mvad)
#' }
build_mmm <-
  function(observations,transition_probs,initial_probs,
           formula, data, coefficients, cluster_names = NULL){

    if (!inherits(observations, "stslist")) {
      stop("The build_mmm function can only be used for single-channel sequence data (as an stslist object). Use the mc_to_sc_data function to convert multiple stslist into single-channel state sequences.")
    }

    n_sequences<-nrow(observations)
    length_of_sequences<-ncol(observations)

    symbol_names <- alphabet(observations)
    n_symbols<-length(symbol_names)

    if (is.list(transition_probs)){
      n_clusters <- length(transition_probs)
    }else{
      stop("Transition_probs is not a list.")
    }
    if(length(initial_probs)!=n_clusters)
      stop("Unequal list lengths of transition_probs and initial_probs.")

    if(is.null(cluster_names)){
      cluster_names <- paste("Cluster", 1:n_clusters)
    }else if(length(cluster_names)!=n_clusters){
      warning("The length of argument cluster_names does not match the number of clusters. Names were not used.")
      cluster_names <- paste("Cluster", 1:n_clusters)
    }

    for(i in 1:n_clusters){
      if (!is.matrix(transition_probs[[i]])) {
        stop(paste("Object provided in transition_probs for cluster", i, "is not a matrix."))
      }
      if (!is.vector(initial_probs[[i]])) {
        stop(paste("Object provided in initial_probs for cluster", i, "is not a vector."))
      }
    }

    # States
    n_states <-rep(n_symbols,n_clusters)

    if (any(rep(n_states, each = 2) != unlist(lapply(transition_probs, dim)))) {
      stop("Transition matrices must be square matrices.")
    }

    state_names <- replicate(n_clusters,symbol_names, simplify = FALSE)


    for(i in 1:n_clusters){
      if (!isTRUE(all.equal(rowSums(transition_probs[[i]]),
                            rep(1, n_states[i]), check.attributes=FALSE))) {
        stop(paste("Row sums of the transition probabilities in cluster", i, "do not sum to one."))
      }
      if (!isTRUE(all.equal(sum(initial_probs[[i]]), 1, check.attributes=FALSE))){
        stop(paste("Initial state probabilities in cluster", i, "do not sum to one."))
      }
    }



    emission_probs <- replicate(n_clusters, diag(n_symbols), simplify = FALSE)
    for(i in 1:n_clusters){
      dimnames(transition_probs[[i]]) <- dimnames(emission_probs[[i]]) <- list(from=state_names[[i]],to=state_names[[i]])
      names(initial_probs[[i]]) <- state_names[[i]]
    }


    if(missing(formula)){
      formula <- stats::formula(rep(1, n_sequences) ~ 1)
    }
    if(missing(data))
      data <- environment(formula)
    if(inherits(formula, "formula")){
      X <- model.matrix(formula, data)
      if(nrow(X)!=n_sequences){
        if(length(all.vars(formula)) > 0 && sum(!complete.cases(data[all.vars(formula)])) > 0){
          stop("Missing cases are not allowed in covariates. Use e.g. the complete.cases function to detect them, then fix, impute, or remove.")
        }else{
          stop("Number of subjects in data for covariates does not match the number of subjects in the sequence data.")
        }
      }
      n_covariates<-ncol(X)
    }else{
      stop("Object given for argument formula is not of class formula.")
    }
    if(missing(coefficients)){
      coefficients<-matrix(0,n_covariates,n_clusters)
    } else {
      if(ncol(coefficients)!=n_clusters | nrow(coefficients)!=n_covariates)
        stop("Wrong dimensions of coefficients.")
      coefficients[,1]<-0
    }


    rownames(coefficients) <- colnames(X)
    colnames(coefficients) <- cluster_names

    names(transition_probs) <- names(initial_probs) <- names(emission_probs) <- cluster_names

      nobs <- sum(!(observations == attr(observations, "nr") |
                      observations == attr(observations, "void") |
                      is.na(observations)))



    model <- structure(list(observations=observations, transition_probs=transition_probs,
                            emission_probs=emission_probs, initial_probs=initial_probs,
                            coefficients=coefficients, X=X, cluster_names=cluster_names, state_names=state_names,
                            symbol_names=symbol_names, channel_names=NULL,
                            length_of_sequences=length_of_sequences,
                            n_sequences=n_sequences, n_clusters=n_clusters,
                            n_symbols=n_symbols, n_states=n_states,
                            n_channels=1,
                            n_covariates=n_covariates, formula = formula), class = "mhmm",
                       nobs = nobs,
                       df = sum(unlist(initial_probs) > 0) - n_clusters + sum(unlist(transition_probs) > 0) - sum(n_states) +
                         n_covariates * (n_clusters - 1),
      type = "mmm")
    model
  }
