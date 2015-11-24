#' Build a Latent Class Model
#' 
#' Function \code{build_lcm} is a shortcut for constructing a latent class model as an restricted case of \code{mhmm} object.
#' 
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param emission_probs A matrix containing emission probabilities for each class by rows, 
#' or in case of multichannel data a list of such matrices. 
#'   Note that the matrices must have dimensions k x s where k is the number of 
#'   latent classes and s is the number of unique symbols (observed states) in the 
#'   data. Emission probabilities should follow the ordering of the alphabet of 
#'   observations (\code{alphabet(observations)}, returned as \code{symbol_names}).
#' @param formula Covariates as an object of class \code{\link{formula}}, 
#' left side omitted.
#' @param data An optional data frame, list or environment containing the variables 
#' in the model. If not found in data, the variables are taken from 
#' \code{environment(formula)}.
#' @param coefficients An optional $k x l$ matrix of regression coefficients for time-constant 
#'   covariates for mixture probabilities, where $l$ is the number of clusters and $k$
#'   is the number of covariates. A logit-link is used for mixture probabilities.
#'   The first column is set to zero.
#' @param cluster_names A vector of optional names for the classes/clusters.
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
#' @seealso \code{\link{fit_mhmm}} for fitting mixture Hidden Markov models; 
#' \code{\link{summary.mhmm}} for a summary of a MHMM; \code{\link{separate_mhmm}} for 
#' reorganizing a MHMM into a list of separate hidden Markov models; and
#' \code{\link{plot.mhmm}} for plotting \code{mhmm} objects.
#' 
#' @examples
#' 
#' data(biofam3c)
#' 
#' ## Building sequence objects
#' marr.seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child.seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left.seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#' 
#' ## Choosing colors
#' attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
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
#'   observations = list(marr.seq, child.seq, left.seq),
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
build_lcm <- 
  function(observations, emission_probs, 
    formula, data, coefficients, cluster_names= NULL, channel_names = NULL){
    
    

    if (is.list(emission_probs)) {
      n_channels <- length(emission_probs)
      n_clusters <- nrow(emission_probs[[1]])
      for (i in 1:n_channels) {
        if (nrow(emission_probs[[i]]) != n_clusters)
          stop("Different number of rows in emission_probs.")
      }
    } else {
      n_channels <- 1
      n_clusters <- nrow(emission_probs)
    }
    
    if(is.null(cluster_names)){
      state_names <- cluster_names <- paste("Class", 1:n_clusters)
    }else if(length(cluster_names)!=n_clusters){
      warning("The length of argument cluster_names does not match the number of clusters. Names were not used.")
      state_names <- cluster_names <- paste("Class", 1:n_clusters)
    }
    
    # States
    n_states <- rep(1, n_clusters)
    
    
    initial_probs <- replicate(n_clusters, 1, simplify = FALSE)
    
    
    
    # Single channel but observations is a list
    if (is.list(observations) && !inherits(observations, "stslist") && length(observations)==1) {
      observations <- observations[[1]]
    }
    
    
    if (n_channels > 1) {
      if (length(observations) != n_channels) {
        stop("Number of channels defined by emission_probs differs from one defined by observations.")
      }
      
      n_sequences <- nrow(observations[[1]])
      length_of_sequences <- ncol(observations[[1]])
      
      emission_probs_list <- vector("list", n_clusters)
      for (i in 1:n_clusters) {
        emission_probs_list[[i]] <- vector("list", n_channels)
        for (j in 1:n_channels) {
          emission_probs_list[[i]][[j]] <- emission_probs[[j]][i, , drop = FALSE]
        }
      }
      emission_probs <- emission_probs_list
      symbol_names <- lapply(observations,alphabet)
      n_symbols <- sapply(symbol_names,length)
      for (i in 1:n_clusters) {
        if (any(n_symbols!=sapply(emission_probs[[i]],ncol))) {
          stop(paste("Number of columns in emission_probs of cluster", i, "is not equal to the number of symbols."))
        }
        if (!isTRUE(all.equal(c(sapply(emission_probs[[i]],rowSums)),
          rep(1,n_channels*n_states[i]),check.attributes=FALSE))) {
          stop(paste("Emission probabilities in emission_probs of cluster", i, "do not sum to one."))
        }
        if (is.null(channel_names)) {
          channel_names<- paste("Channel", 1:n_channels)
        } else if (length(channel_names)!=n_channels) {
          warning("The length of argument channel_names does not match the number of channels. Names were not used.")
          channel_names<- paste("Channel", 1:n_channels)
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
    
    transition_probs <- replicate(n_clusters, matrix(1), simplify = FALSE)
    
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
      df = sum(unlist(initial_probs) > 0) - n_clusters + 
        sum(unlist(emission_probs) > 0) - sum(n_states) * n_channels + n_covariates * (n_clusters - 1),
      type = "lcm")
    model
  }
