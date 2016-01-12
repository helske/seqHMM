#' Build a Latent Class Model
#' 
#' Function \code{build_lcm} is a shortcut for constructing a latent class model 
#' as a restricted case of an \code{mhmm} object.
#' 
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param emission_probs A matrix containing emission probabilities for each class by rows, 
#'   or in case of multichannel data a list of such matrices. 
#'   Note that the matrices must have dimensions k x s where k is the number of 
#'   latent classes and s is the number of unique symbols (observed states) in the 
#'   data. Emission probabilities should follow the ordering of the alphabet of 
#'   observations (\code{alphabet(observations)}, returned as \code{symbol_names}).
#' @param formula Covariates as an object of class \code{\link{formula}}, 
#' left side omitted.
#' @param data An optional data frame, list or environment containing the variables 
#' in the model. If not found in data, the variables are taken from 
#' \code{environment(formula)}.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for 
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number 
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the classes/clusters.
#' @param channel_names A vector of optional names for the channels.
#' @return Object of class \code{mhmm} with the following elements:
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
#' # Simulate observations from two classes
#' set.seed(123)
#' obs <- seqdef(rbind(
#'   matrix(sample(letters[1:3], 5000, TRUE, prob = c(0.1, 0.6, 0.3)), 500, 10),
#'   matrix(sample(letters[1:3], 2000, TRUE, prob = c(0.4, 0.4, 0.2)), 200, 10)))
#' 
#' # Initialize the model
#' model <- build_lcm(obs, 
#'   # Simulating emission probs (matrix of n_clusters rows and n_symbols columns)
#'   # Note that the simulation function is made for HMMs -> number of clusters defined with n_states
#'   simulate_emission_probs(n_states = 2, n_symbols = 3))
#' 
#' # Estimate model parameters
#' fit <- fit_model(model)
#' 
#' # How many of the observations were correctly classified:
#' sum(summary(fit$model)$most_probable_cluster == rep(c("Class 2", "Class 1"), times = c(500, 200)))
#' 
#' ############################################################
#' 
#' # LCM for longitudinal data
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
#' # Starting values for emission probabilities
#' emiss <- rbind(rep(1/6, 6), rep(1/6, 6))
#' 
#' # Initialize the LCM
#' init_lcm_mvad <- build_lcm(observations = mvad_seq,
#'   emission_probs = emiss, formula = ~male, data = mvad)
#'
#' # Estimate model parameters (EM algorithm with random restarts)
#' lcm_mvad <- fit_model(init_lcm_mvad, 
#'   control_em = list(restart = list(times = 5)))$model
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
#' model <- build_lcm(seqdef(birthwt$low), diag(2), ~ age + lwt + smoke + ht, birthwt)
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
#' X <- cbind(1, x1 = runif(n, 0, 1), x2 =  runif(n, 0, 1))
#' coefs <- cbind(0,c(-2, 5, -2), c(0, -2, 2))
#' pr <- exp(X %*% coefs)  + rnorm(n*3)
#' pr <- pr/rowSums(pr)
#' y <- apply(pr, 1, which.max)
#' table(y)
#' 
#' model <- build_lcm(seqdef(y), diag(3), ~ x1 + x2,  data.frame(X[, -1]))
#' fit <- fit_model(model)
#' summary(fit$model)
#' summary(multinom(y ~ x1 + x2, data = data.frame(X[,-1])))
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
        stop("Number of channels defined by emission_probs differs from the one defined by observations. Give emission_probs as a matrix with n_cluster rows and n_symbol columns (for single-channel data) or a list of such matrices (for multichannel data).")
      }
      
      if (length(unique(sapply(observations, nrow))) > 1) {
        stop("The number of subjects (rows) is not the same in all channels.")
      }
      if (length(unique(sapply(observations, ncol))) > 1) {
        stop("The length of the sequences (number of columns) is not the same in all channels.")
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
      n_symbols <- lengths(symbol_names)
      for (i in 1:n_clusters) {
        if (any(n_symbols!=sapply(emission_probs[[i]],ncol))) {
          stop(paste("Number of columns in emission_probs of cluster", i, "is not equal to the number of symbols."))
        }
        if (!isTRUE(all.equal(c(sapply(emission_probs[[i]],rowSums)),
          rep(1,n_channels*n_states[i]),check.attributes=FALSE))) {
          stop(paste("Emission probabilities in emission_probs of cluster", i, "do not sum to one."))
        }
        if (is.null(channel_names)) {
          if(is.null(channel_names <- names(observations))){
            channel_names <- paste("Channel", 1:n_channels)
          }
        } else if (length(channel_names)!=n_channels) {
          warning("The length of argument channel_names does not match the number of channels. Names were not used.")
          channel_names<- paste("Channel", 1:n_channels)
        }
        for (j in 1:n_channels) {
          dimnames(emission_probs[[i]][[j]])<-list(state_names=state_names[[i]],symbol_names=symbol_names[[j]])
        }
        names(emission_probs[[i]])<-channel_names
        
      }
    } else {
      
      if (is.null(channel_names)) {
        channel_names <- "Observations"
      }      
      n_sequences<-nrow(observations)
      length_of_sequences<-ncol(observations)
      symbol_names<-alphabet(observations)
      n_symbols<-length(symbol_names)
      
      emission_probs_list <- vector("list", n_clusters)
      for (i in 1:n_clusters) {
        emission_probs_list[[i]] <- emission_probs[i, , drop = FALSE]
      }
      emission_probs <- emission_probs_list
      
      for(i in 1:n_clusters){
        if(n_symbols!=ncol(emission_probs[[i]]))
          stop(paste("Number of columns in emission_probs of cluster", i, "is not equal to the number of symbols."))
        if(!isTRUE(all.equal(rep(1,n_states[i]),rowSums(emission_probs[[i]]),check.attributes=FALSE)))
          stop(paste("Emission probabilities in emission_probs of cluster", i, "do not sum to one."))
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
    
    for (i in 1:n_clusters) {
      names(initial_probs[[i]]) <- rownames(transition_probs[[i]]) <- 
        colnames(transition_probs[[i]]) <- state_names[[i]]
    }
    
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
      df = sum(unlist(emission_probs) > 0) - sum(n_states) * n_channels + 
        n_covariates * (n_clusters - 1),
      type = "lcm")
    model
  }
