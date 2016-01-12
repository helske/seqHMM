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
#' @param state_names A list of optional labels for the hidden states. If \code{NULL},
#' the state names are taken as row names of transition matrices. If this is also \code{NULL},
#' numbered states are used.
#' @param channel_names A vector of optional names for the channels.
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
#' @seealso \code{\link{fit_model}} for fitting mixture Hidden Markov models;
#' \code{\link{summary.mhmm}} for a summary of a MHMM; \code{\link{separate_mhmm}} for
#' reorganizing a MHMM into a list of separate hidden Markov models; and
#' \code{\link{plot.mhmm}} for plotting \code{mhmm} objects.
#'
#' @examples
#'
#' data("biofam3c")
#'
#' ## Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' ## Choosing colors
#' attr(marr_seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(child_seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(left_seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#'
#' ## Starting values for emission probabilities
#'
#' # Cluster 1
#' B1_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.3, 0.6, 0.1, # High probability for married
#'     0.3, 0.3, 0.4), # High probability for divorced
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B1_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.9, 0.1),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B1_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9, # High probability for having left home
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 2
#'
#' B2_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.7, 0.2, 0.1),
#'   nrow = 4, ncol = 3, byrow = TRUE)
#'
#' B2_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.9, 0.1,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' B2_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 4, ncol = 2, byrow = TRUE)
#'
#' # Cluster 3
#' B3_marr <- matrix(
#'   c(0.8, 0.1, 0.1, # High probability for single
#'     0.8, 0.1, 0.1,
#'     0.8, 0.1, 0.1,
#'     0.1, 0.8, 0.1, # High probability for married
#'     0.3, 0.4, 0.3,
#'     0.1, 0.1, 0.8), # High probability for divorced
#'   nrow = 6, ncol = 3, byrow = TRUE)
#'
#' B3_child <- matrix(
#'   c(0.9, 0.1, # High probability for childless
#'     0.9, 0.1,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#'
#' B3_left <- matrix(
#'   c(0.9, 0.1, # High probability for living with parents
#'     0.1, 0.9,
#'     0.5, 0.5,
#'     0.5, 0.5,
#'     0.1, 0.9,
#'     0.1, 0.9),
#'   nrow = 6, ncol = 2, byrow = TRUE)
#'
#' # Starting values for transition matrices
#' A1 <- matrix(
#'   c(0.80, 0.16, 0.03, 0.01,
#'     0,    0.90, 0.07, 0.03,
#'     0,    0,    0.90, 0.10,
#'     0,    0,    0,       1),
#'   nrow = 4, ncol = 4, byrow = TRUE)
#'
#' A2 <- matrix(
#'   c(0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
#'     0,    0.70, 0.10, 0.10, 0.05, 0.05,
#'     0,    0,    0.85, 0.01, 0.10, 0.04,
#'     0,    0,    0,    0.90, 0.05, 0.05,
#'     0,    0,    0,    0,    0.90, 0.10,
#'     0,    0,    0,    0,    0,       1),
#'   nrow = 6, ncol = 6, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#'
#' # Birth cohort
#' biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr, c(1908, 1935, 1945, 1957))
#' biofam3c$covariates$cohort <- factor(
#'   biofam3c$covariates$cohort, labels=c("1909-1935", "1936-1945", "1946-1957"))
#'
#' # Build mixture HMM
#' init_mhmm_bf <- build_mhmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   transition_probs = list(A1, A1, A2),
#'   emission_probs = list(list(B1_marr, B1_child, B1_left),
#'     list(B2_marr, B2_child, B2_left),
#'     list(B3_marr, B3_child, B3_left)),
#'   formula = ~sex + cohort, data = biofam3c$covariates,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Marriage", "Parenthood", "Residence"),
#'   state_names = list(paste("State", 1:4), paste("State", 1:4),
#'                      paste("State", 1:6)))
#'
build_mhmm <-
  function(observations,transition_probs,emission_probs,initial_probs,
           formula, data, coefficients, cluster_names=NULL, state_names=NULL, channel_names=NULL){

    if (is.list(transition_probs)){
      n_clusters<-length(transition_probs)
    }else{
      stop("Transition_probs is not a list.")
    }
    if(length(emission_probs)!=n_clusters || length(initial_probs)!=n_clusters)
      stop("Unequal list lengths of transition_probs, emission_probs and initial_probs.")

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
    n_states <- unlist(lapply(transition_probs,nrow))

    if (any(rep(n_states, each = 2) != unlist(lapply(transition_probs, dim)))) {
      stop("Transition matrices must be square matrices.")
    }

    if (is.null(state_names)) {
      state_names <- vector("list", n_clusters)
      for(m in 1:n_clusters){
        if (is.null(state_names[[m]] <- rownames(transition_probs[[m]]))) {
          state_names[[m]] <-  paste("State", 1:n_states[m])
        }
      }
    } else {
      for (m in 1:n_clusters) {
        if (length(state_names[[m]]) != n_states[m]) {
          stop(paste0("Length of state_names for cluster ", m, " is not equal to the number of hidden states."))
        }
      }
    }

    for(i in 1:n_clusters){
      if (!isTRUE(all.equal(rowSums(transition_probs[[i]]),
                            rep(1, n_states[i]), check.attributes=FALSE))) {
        stop(paste("Row sums of the transition probabilities in cluster", i, "do not sum to one."))
      }
      if (!isTRUE(all.equal(sum(initial_probs[[i]]), 1, check.attributes=FALSE))){
        stop(paste("Initial state probabilities in cluster", i, "do not sum to one."))
      }
    }




    for(i in 1:n_clusters){
      dimnames(transition_probs[[i]]) <- list(from=state_names[[i]],to=state_names[[i]])
      # Single channel but emission_probs is list of lists
      if(is.list(emission_probs[[i]]) && length(emission_probs[[i]])==1)
        emission_probs[[i]] <- emission_probs[[i]][[1]]
    }




    # Single channel but observations is a list
    if (is.list(observations) && !inherits(observations, "stslist") && length(observations)==1) {
      observations <- observations[[1]]
    }

    n_channels <- ifelse(is.list(emission_probs[[1]]), length(emission_probs[[1]]), 1)

    for(i in 1:n_clusters){
      if (n_channels == 1){
        if (!is.matrix(emission_probs[[i]])) {
          stop(paste("Object provided in emission_probs for cluster", i, "is not a matrix."))
        }
      } else {
        for (j in 1:n_channels){
          if (!is.matrix(emission_probs[[i]][[j]])) {
            stop(paste("Object provided in emission_probs for cluster", i, "and channel", j, "is not a matrix."))
          }
        }
      }
    }

    if (n_channels>1 && any(sapply(emission_probs,length)!=n_channels)) {
      stop("Number of channels defined by emission matrices differ from each other.")
    }

    if(n_channels>1){
      if(length(observations)!=n_channels){
        stop("Number of channels defined by emission_probs differs from one defined by observations.")
      }

      if (length(unique(sapply(observations, nrow))) > 1) {
        stop("The number of subjects (rows) is not the same in all channels.")
      }
      if (length(unique(sapply(observations, ncol))) > 1) {
        stop("The length of the sequences (number of columns) is not the same in all channels.")
      }

      n_sequences <- nrow(observations[[1]])
      length_of_sequences <- ncol(observations[[1]])


      symbol_names <- lapply(observations, alphabet)
      n_symbols <- lengths(symbol_names)
      for (i in 1:n_clusters) {
        if (length(initial_probs[[i]]) != n_states[i]) {
          stop(paste("Length of initial_probs of cluster", i, "is not equal to the number of states."))
        }
        if (any(lapply(emission_probs[[i]],nrow) != n_states[i])) {
          stop(paste("Number of rows in emission_probs of cluster", i, "is not equal to the number of states."))
        }

        if (any(n_symbols != sapply(emission_probs[[i]],ncol))) {
          stop(paste("Number of columns in emission_probs of cluster", i, "is not equal to the number of symbols."))
        }
        if (!isTRUE(all.equal(c(sapply(emission_probs[[i]], rowSums)),
                              rep(1, n_channels * n_states[i]), check.attributes = FALSE))) {
          stop(paste("Emission probabilities in emission_probs of cluster", i, "do not sum to one."))
        }
        if (is.null(channel_names)) {
          if(is.null(channel_names <- names(observations))){
            channel_names <- paste("Channel", 1:n_channels)
          }
        } else if (length(channel_names)!=n_channels) {
          warning("The length of argument channel_names does not match the number of channels. Names were not used.")
          channel_names<-paste("Channel", 1:n_channels)
        }
        for (j in 1:n_channels) {
          dimnames(emission_probs[[i]][[j]])<-list(state_names=state_names[[i]],symbol_names=symbol_names[[j]])
        }
        names(emission_probs[[i]])<-channel_names
        names(initial_probs[[i]]) <- state_names[[i]]
      }
    } else {
      n_channels <- 1
      if (is.null(channel_names)) {
        channel_names <- "Observations"
      }
      n_sequences<-nrow(observations)
      length_of_sequences<-ncol(observations)
      symbol_names<-alphabet(observations)
      n_symbols<-length(symbol_names)

      for(i in 1:n_clusters){
        if(n_states[i]!=dim(emission_probs[[i]])[1])
          stop("Number of rows in emission_probs is not equal to the number of states.")
        if(n_symbols!=dim(emission_probs[[i]])[2])
          stop("Number of columns in emission_probs is not equal to the number of symbols.")
        if(!isTRUE(all.equal(rep(1,n_states[i]),rowSums(emission_probs[[i]]),check.attributes=FALSE)))
          stop("Emission probabilities in emission_probs do not sum to one.")
        dimnames(emission_probs[[i]])<-list(state_names=state_names[[i]],symbol_names=symbol_names)
        names(initial_probs[[i]]) <- state_names[[i]]
      }

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

    names(transition_probs) <- names(emission_probs) <- names(initial_probs) <- cluster_names
    if(n_channels > 1){
      nobs <- sum(sapply(observations, function(x) sum(!(x == attr(observations[[1]], "nr") |
                                                           x == attr(observations[[1]], "void") |
                                                           is.na(x)))))/n_channels
    } else {
      nobs <- sum(!(observations == attr(observations, "nr") |
                      observations == attr(observations, "void") |
                      is.na(observations)))
    }
    model <- structure(list(observations=observations, transition_probs=transition_probs,
                            emission_probs=emission_probs, initial_probs=initial_probs,
                            coefficients=coefficients, X=X, cluster_names=cluster_names, state_names=state_names,
                            symbol_names=symbol_names, channel_names=channel_names,
                            length_of_sequences=length_of_sequences,
                            n_sequences=n_sequences, n_clusters=n_clusters,
                            n_symbols=n_symbols, n_states=n_states,
                            n_channels=n_channels,
                            n_covariates=n_covariates, formula = formula), class = "mhmm",
                       nobs = nobs,
                       df = sum(unlist(initial_probs) > 0) - n_clusters + sum(unlist(transition_probs) > 0) - sum(n_states) +
                         sum(unlist(emission_probs) > 0) - sum(n_states) * n_channels + n_covariates * (n_clusters - 1),
      type = "mhmm")
    model
  }
