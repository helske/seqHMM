#' Build a Mixture Hidden Markov Model
#' 
#' Function build_mhmm constructs a mixture of hidden Markov models.
#' 
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param transition_matrix A list of matrices of transition 
#'   probabilities for submodels of each cluster.
#' @param emission_matrix A list which contains matrices of emission probabilities or
#'   a list of such objects (one for each channel) for submodels of each cluster. 
#'   Note that the matrices must have dimensions m x s where m is the number of 
#'   hidden states and s is the number of unique symbols (observed states) in the 
#'   data.
#' @param initial_probs A list which contains vectors of initial state 
#'   probabilities for submodels of each cluster.
#' @param formula Covariates as an object of class \code{\link{formula}}, 
#' left side omitted.
#' @param data An optional data frame, list or environment containing the variables 
#' in the model. If not found in data, the variables are taken from 
#' \code{environment(formula)}.
#' @param beta An optional k x l matrix of regression coefficients for time-constant 
#'   covariates for mixture probabilities, where l is the number of clusters and k
#'   is the number of covariates. A logit-link is used for mixture probabilities.
#'   The first column is set to zero.
#' @param cluster_names A vector of optional names for the clusters.
#' @param state_names A list of optional labels for the hidden states.
#' @param channel_names A vector of optional names for the channels.
#' @return Object of class \code{mhmm}.
#' @seealso \code{\link{fit_mhmm}} for fitting mixture Hidden Markov models.
#' 
#' @examples
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[complete.cases(biofam[c(2:4)]),]
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 5) > 0) | 
#'             (rowSums(bf == 7) > 0 & rowSums(bf == 6) > 0),]
#' children[rownames(bf) %in% rownames(div) & bf == 7] <- "Children"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' married[bf == 7] <- "Divorced"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 & rowSums(bf == 3) == 0 &  
#'           rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0) | 
#'          (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 & rowSums(bf == 3) == 0 &  
#'          rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0),]
#' left[rownames(bf) %in% rownames(wp) & bf == 7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start = 15)
#' marr.seq <- seqdef(married, start = 15)
#' left.seq <- seqdef(left, start = 15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' B2_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.01, 0.99), nrow = 4, ncol = 2, byrow = TRUE)
#'                      
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.29, 0.7, 0.01),
#'                    nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE) 
#' 
#' # Cluster 3
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.01, 0.99), nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow = 6, ncol = 3, byrow = TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.50, 0.50,
#'                     0.01, 0.99,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 6, ncol = 2, byrow = TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                  0,    0.9, 0.07, 0.03, 
#'                  0,      0,  0.9,  0.1, 
#'                  0,      0,    0,    1), 
#'              nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                  0,  0.7,  0.1,   0.1, 0.05, 0.05,
#'                  0,    0, 0.85,  0.01,  0.1, 0.04,
#'                  0,    0,    0,   0.9, 0.05, 0.05,
#'                  0,    0,    0,     0,  0.9,  0.1,
#'                  0,    0,    0,     0,    0,    1), 
#'              nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Creating covariate swiss
#' biofam$swiss <- biofam$nat_1_02 == "Switzerland"
#' biofam$swiss[biofam$swiss == TRUE] <- "Swiss"
#' biofam$swiss[biofam$swiss == FALSE] <- "Other"
#' 
#' # Birth cohort
#' biofam$cohort <- cut(biofam$birthyr, c(1908, 1935, 1945, 1957))
#' biofam$cohort <- factor(
#'   biofam$cohort, labels=c("1909-1935", "1936-1945", "1946-1957")
#' )
#' 
#' # Build mixture HMM
#' bMHMM <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = list(A1,A1,A2),
#'   emission_matrix = list(list(B1_child, B1_marr, B1_left),
#'                         list(B2_child, B2_marr, B2_left), 
#'                         list(B3_child, B3_marr, B3_left)),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   formula = ~ sex + cohort, data = biofam,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Left home")
#'   )
#'                     
build_mhmm <- 
  function(observations,transition_matrix,emission_matrix,initial_probs, 
    formula, data, beta, cluster_names=NULL, state_names=NULL, channel_names=NULL){
    
    n_clusters<-length(transition_matrix)
    if(length(emission_matrix)!=n_clusters || length(initial_probs)!=n_clusters)
      stop("Unequal lengths of transition_matrix, emission_matrix and initial_probs.")
    
    if(is.null(cluster_names)){
      cluster_names <- paste("Cluster", 1:n_clusters)
    }else if(length(cluster_names)!=n_clusters){
      warning("The length of argument cluster_names does not match the number of clusters. Names were not used.")
      cluster_names <- paste("Cluster", 1:n_clusters)
    }
    
    model <- vector("list", length = n_clusters)
    
    # States
    n_states <- unlist(lapply(transition_matrix,nrow))
    
    if(any(rep(n_states,each=2)!=unlist(lapply(transition_matrix,dim))))
      stop("Transition matrices must be square matrices.")
    
    if(is.null(state_names)){
      state_names <- vector("list", n_clusters)
      for(m in 1:n_clusters){
        state_names[[m]] <- as.character(1:n_states[m])
      }
    }
    
    if(!all(1==unlist(sapply(transition_matrix,rowSums))))
      stop("Transition probabilities in transition_matrix do not sum to one.")
    
    if(!all(1==unlist(sapply(initial_probs,sum))))
      stop("Initial state probabilities do not sum to one.")
    
    for(i in 1:n_clusters){
      
      dimnames(transition_matrix[[i]]) <- list(from=state_names[[i]],to=state_names[[i]])
      # Single channel but emission_matrix is list of lists  
      if(is.list(emission_matrix[[i]]) && length(emission_matrix[[i]])==1)   
        emission_matrix[[i]] <- emission_matrix[[i]][[1]]
    }
    
    
    
    
    # Single channel but observations is a list
    if(is.list(observations) && !inherits(observations, "stslist") && length(observations)==1)
      observations <- observations[[1]]
    
    n_channels <- ifelse(is.list(emission_matrix[[1]]),length(emission_matrix[[1]]),1)
    
    if(n_channels>1 && any(sapply(emission_matrix,length)!=n_channels))
      stop("Number of channels defined by emission matrices differ from each other.")
    
    if(n_channels>1){
      if(length(observations)!=n_channels){
        stop("Number of channels defined by emission_matrix differs from one defined by observations.")
      }
      
      
      n_sequences<-nrow(observations[[1]])
      length_of_sequences<-ncol(observations[[1]])
      
      
      symbol_names<-lapply(observations,alphabet)
      n_symbols<-sapply(symbol_names,length)
      for(i in 1:n_clusters){
        if(any(lapply(emission_matrix[[i]],nrow)!=n_states[i]))
          stop(paste("Number of rows in emission_matrix of cluster", i, "is not equal to the number of states."))
        
        if(any(n_symbols!=sapply(emission_matrix[[i]],ncol)))
          stop(paste("Number of columns in emission_matrix of cluster", i, "is not equal to the number of symbols."))
        if(!isTRUE(all.equal(c(sapply(emission_matrix[[i]],rowSums)),
          rep(1,n_channels*n_states[i]),check.attributes=FALSE)))
          stop(paste("Emission probabilities in emission_matrix of cluster", i, "do not sum to one."))
        if(is.null(channel_names)){
          channel_names<-as.character(1:n_channels)
        }else if(length(channel_names)!=n_channels){
          warning("The length of argument channel_names does not match the number of channels. Names were not used.")
          channel_names<-as.character(1:n_channels)
        }
        for(j in 1:n_channels)
          dimnames(emission_matrix[[i]][[j]])<-list(state_names=state_names[[i]],symbol_names=symbol_names[[j]])
        names(emission_matrix[[i]])<-channel_names
        names(initial_probs[[i]]) <- state_names[[i]]
      }
    } else {
      n_channels <- 1
      channel_names<-NULL
      n_sequences<-nrow(observations)
      length_of_sequences<-ncol(observations)
      symbol_names<-alphabet(observations)
      n_symbols<-length(symbol_names)
      
      for(i in 1:n_clusters){
        if(n_states[i]!=dim(emission_matrix[[i]])[1])
          stop("Number of rows in emission_matrix is not equal to the number of states.")
        if(n_symbols!=dim(emission_matrix[[i]])[2])
          stop("Number of columns in emission_matrix is not equal to the number of symbols.")
        if(!isTRUE(all.equal(rep(1,n_states[i]),rowSums(emission_matrix[[i]]),check.attributes=FALSE)))
          stop("Emission probabilities in emission_matrix do not sum to one.")
        dimnames(emission_matrix[[i]])<-list(state_names=state_names[[i]],symbol_names=symbol_names)
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
    if(missing(beta)){
      beta<-matrix(0,n_covariates,n_clusters)
    } else {
      if(ncol(beta)!=n_clusters | nrow(beta)!=n_covariates)
        stop("Wrong dimensions of beta.")
      beta[,1]<-0
    }       
    
    
    rownames(beta) <- colnames(X)
    colnames(beta) <- cluster_names
    
    names(transition_matrix) <- names(emission_matrix) <- names(initial_probs) <- cluster_names
    if(n_channels > 1){
    nobs <- sum(sapply(observations, function(x) sum(!(x == attr(observations[[1]], "nr") |
        x == attr(observations[[1]], "void") |
        is.na(x)))))/n_channels
    } else {
      nobs <- sum(!(x == attr(observations, "nr") |
          x == attr(observations, "void") |
          is.na(x)))
    }
    model <- structure(list(observations=observations, transition_matrix=transition_matrix,
      emission_matrix=emission_matrix, initial_probs=initial_probs,
      beta=beta, X=X, cluster_names=cluster_names, state_names=state_names, 
      symbol_names=symbol_names, channel_names=channel_names, 
      length_of_sequences=length_of_sequences,
      n_sequences=n_sequences, n_clusters=n_clusters,
      n_symbols=n_symbols, n_states=n_states,
      n_channels=n_channels,
      n_covariates=n_covariates), class = "mhmm", 
      nobs = nobs,
      df = sum(unlist(initial_probs) > 0) - n_clusters + sum(unlist(transition_matrix) > 0) - sum(n_states) + 
        sum(unlist(emission_matrix) > 0) - sum(n_states) * n_channels + n_covariates * (n_clusters - 1))
    model
  }