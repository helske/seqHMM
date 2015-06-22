#' Build a Mixture Hidden Markov Model
#' 
#' Function buildMixHMM constructs a mixture of hidden Markov models.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param transitionMatrix A list of matrices of transition 
#'   probabilities for submodels of each cluster.
#' @param emissionMatrix A list which contains matrices of emission probabilities or
#'   a list of such objects (one for each channel) for submodels of each cluster. 
#'   Note that the matrices must have dimensions m x s where m is the number of 
#'   hidden states and s is the number of unique symbols (observed states) in the 
#'   data.
#' @param initialProbs A list which contains vectors of initial state 
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
#' @param clusterNames A vector of optional names for the clusters.
#' @param stateNames A list of optional labels for the hidden states.
#' @param channelNames A vector of optional names for the channels.
#' @return Object of class \code{mixHMModel}
#' @seealso \code{\link{fitMixHMM}} for fitting mixture Hidden Markov models.
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
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6 | bf==7
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf==7)>0 & rowSums(bf==5)>0) | 
#'             (rowSums(bf==7)>0 & rowSums(bf==6)>0),]
#' children[rownames(bf) %in% rownames(div) & bf==7] <- "Children"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' married[bf==7] <- "Divorced"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf==7)>0 & rowSums(bf==2)>0 & rowSums(bf==3)==0 &  rowSums(bf==5)==0 &  rowSums(bf==6)==0) | 
#'            (rowSums(bf==7)>0 & rowSums(bf==4)>0 & rowSums(bf==3)==0 &  rowSums(bf==5)==0 &  rowSums(bf==6)==0),]
#' left[rownames(bf) %in% rownames(wp) & bf==7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start=15)
#' marr.seq <- seqdef(married, start=15)
#' left.seq <- seqdef(left, start=15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01
#'                     0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.98, 0.01, # High probability for married
#'                      0.29, 0.7, 0.01),
#'                    nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=4, ncol=2, byrow=TRUE) 
#' 
#' # Sinkkuvanhemmat ja kotona asuvat yhdessä
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                       0.99, 0.01,
#'                       0.01, 0.99,
#'                       0.99, 0.01,
#'                       0.01, 0.99,
#'                       0.01, 0.99), nrow=6, ncol=2, byrow=TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.98, 0.01,
#'                      0.01, 0.98, 0.01, # High probability for married
#'                      0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow=6, ncol=3, byrow=TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                      0.99, 0.01,
#'                      0.50, 0.50,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=6, ncol=2, byrow=TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                0,    0.9, 0.07, 0.03, 
#'                0,      0,  0.9,  0.1, 
#'                0,      0,    0,    1), 
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                0,    0.7,  0.1,   0.1, 0.05, 0.05,
#'                0,      0,  0.85, 0.01,  0.1, 0.04,
#'                0,      0,    0,   0.9, 0.05, 0.05,
#'                0,      0,    0,     0,  0.9,  0.1,
#'                0,      0,    0,     0,    0,    1), 
#'              nrow=6, ncol=6, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initialProbs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initialProbs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Creating covariate swiss
#' bio$swiss <- bio$nat_1_02=="Switzerland"
#' bio$swiss[bio$swiss==TRUE] <- "Swiss"
#' bio$swiss[bio$swiss==FALSE] <- "Other"
#' 
#' # Build mixture HMM
#' bmHMM <- buildMixHMM(observations=list(child.seq, marr.seq, left.seq), 
#'                        transitionMatrix=list(A1,A2,A1), 
#'                        emissionMatrix=list(list(B1_child, B1_marr, B1_left),
#'                                            list(B2_child, B2_marr, B2_left),
#'                                            list(B3_child, B3_marr, B3_left)),
#'                        initialProbs=list(initialProbs1, initialProbs2,
#'                                          initialProbs1), 
#'                        formula=~sex*birthyr+sex*swiss, data=bio,
#'                        clusterNames=c("Cluster 1", "Cluster 2", "Cluster 3"),
#'                        channelNames=c("Parenthood", "Marriage", "Left home"))
#'                     
buildMixHMM <- 
  function(observations,transitionMatrix,emissionMatrix,initialProbs, 
           formula, data, beta, clusterNames=NULL, stateNames=NULL, channelNames=NULL){
    
    numberOfClusters<-length(transitionMatrix)
    if(length(emissionMatrix)!=numberOfClusters || length(initialProbs)!=numberOfClusters)
      stop("Unequal lengths of transitionMatrix, emissionMatrix and initialProbs.")
    
    if(is.null(clusterNames)){
      clusterNames <- paste("Cluster", 1:numberOfClusters)
    }else if(length(clusterNames)!=numberOfClusters){
      warning("The length of argument clusterNames does not match the number of clusters. Names were not used.")
      clusterNames <- paste("Cluster", 1:numberOfClusters)
    }
      
    model <- vector("list", length = numberOfClusters)
    
    # States
    numberOfStates <- unlist(lapply(transitionMatrix,nrow))
    
    if(any(rep(numberOfStates,each=2)!=unlist(lapply(transitionMatrix,dim))))
      stop("Transition matrices must be square matrices.")
    
    if(is.null(stateNames)){
      stateNames <- vector("list", numberOfClusters)
      for(m in 1:numberOfClusters){
        stateNames[[m]] <- as.character(1:numberOfStates[m])
      }
    }
    
    if(!all(1==unlist(sapply(transitionMatrix,rowSums))))
      stop("Transition probabilities in transitionMatrix do not sum to one.")
    
    if(!all(1==unlist(sapply(initialProbs,sum))))
      stop("Initial state probabilities do not sum to one.")

    for(i in 1:numberOfClusters){

      dimnames(transitionMatrix[[i]]) <- list(from=stateNames[[i]],to=stateNames[[i]])
      # Single channel but emissionMatrix is list of lists  
      if(is.list(emissionMatrix[[i]]) && length(emissionMatrix[[i]])==1)   
        emissionMatrix[[i]] <- emissionMatrix[[i]][[1]]
    }
    
    
    
    
    # Single channel but observations is a list
    if(is.list(observations) && length(observations)==1)
      observations <- observations[[1]]
    
    numberOfChannels <- ifelse(is.list(emissionMatrix[[1]]),length(emissionMatrix[[1]]),1)
    
    if(numberOfChannels>1 && any(sapply(emissionMatrix,length)!=numberOfChannels))
      stop("Number of channels defined by emission matrices differ from each other.")
    
    if(numberOfChannels>1){
      if(length(observations)!=numberOfChannels){
        stop("Number of channels defined by emissionMatrix differs from one defined by observations.")
      }
      
      
      numberOfSequences<-nrow(observations[[1]])
      lengthOfSequences<-ncol(observations[[1]])
      

      symbolNames<-lapply(observations,alphabet)
      numberOfSymbols<-sapply(symbolNames,length)
      for(i in 1:numberOfClusters){
        if(any(lapply(emissionMatrix[[i]],nrow)!=numberOfStates[i]))
          stop(paste("Number of rows in emissionMatrix of cluster", i, "is not equal to the number of states."))
        
        if(any(numberOfSymbols!=sapply(emissionMatrix[[i]],ncol)))
          stop(paste("Number of columns in emissionMatrix of cluster", i, "is not equal to the number of symbols."))
        if(!isTRUE(all.equal(c(sapply(emissionMatrix[[i]],rowSums)),
                             rep(1,numberOfChannels*numberOfStates[i]),check.attributes=FALSE)))
          stop(paste("Emission probabilities in emissionMatrix of cluster", i, "do not sum to one."))
        if(is.null(channelNames)){
          channelNames<-as.character(1:numberOfChannels)
        }else if(length(channelNames)!=numberOfChannels){
          warning("The length of argument channelNames does not match the number of channels. Names were not used.")
          channelNames<-as.character(1:numberOfChannels)
        }
        for(j in 1:numberOfChannels)
          dimnames(emissionMatrix[[i]][[j]])<-list(stateNames=stateNames[[i]],symbolNames=symbolNames[[j]])
        names(emissionMatrix[[i]])<-channelNames
      }
    } else {
      numberOfChannels <- 1
      channelNames<-NULL
      numberOfSequences<-nrow(observations)
      lengthOfSequences<-ncol(observations)
      symbolNames<-alphabet(observations)
      numberOfSymbols<-length(symbolNames)
      
      for(i in 1:numberOfClusters){
        if(numberOfStates[i]!=dim(emissionMatrix[[i]])[1])
          stop("Number of rows in emissionMatrix is not equal to the number of states.")
        if(numberOfSymbols!=dim(emissionMatrix[[i]])[2])
          stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
        if(!isTRUE(all.equal(rep(1,numberOfStates[i]),rowSums(emissionMatrix[[i]]),check.attributes=FALSE)))
          stop("Emission probabilities in emissionMatrix do not sum to one.")
        dimnames(emissionMatrix[[i]])<-list(stateNames=stateNames[[i]],symbolNames=symbolNames)
      }
      
    }
    
    
    if(!missing(formula)){
      if(inherits(formula, "formula")){
      X <- model.matrix(formula, data) #[,-1,drop=FALSE]
      if(nrow(X)!=numberOfSequences)
        stop("Number of subjects in data for covariates does not match the number of subjects in the sequence data.")
      numberOfCovariates<-ncol(X)
      }else{
        stop("Object given for argument formula is not of class formula.")
      }
      if(missing(beta)){
        beta<-matrix(0,numberOfCovariates,numberOfClusters)
      } else {
        if(ncol(beta)!=numberOfClusters | nrow(beta)!=numberOfCovariates)
          stop("Wrong dimensions of beta.")
        beta[,1]<-0
      }       
    } else { #Just intercept
      numberOfCovariates <-1
      X <- matrix(1,nrow=numberOfSequences)
      beta <- matrix(0,1,numberOfClusters)        
    }
    
    rownames(beta) <- colnames(X)
    colnames(beta) <- clusterNames
    
    names(transitionMatrix) <- names(emissionMatrix) <- names(initialProbs) <- clusterNames
    
    pr <- exp(X%*%beta)
    clusterProbabilities <- pr/rowSums(pr)
    
    model<-list(observations=observations, transitionMatrix=transitionMatrix,
                emissionMatrix=emissionMatrix, initialProbs=initialProbs,
                beta=beta, X=X, clusterNames=clusterNames, stateNames=stateNames, 
                symbolNames=symbolNames, channelNames=channelNames, 
                lengthOfSequences=lengthOfSequences,
                numberOfSequences=numberOfSequences, numberOfClusters=numberOfClusters,
                numberOfSymbols=numberOfSymbols, numberOfStates=numberOfStates,
                numberOfChannels=numberOfChannels,
                numberOfCovariates=numberOfCovariates, 
                clusterProbabilities=clusterProbabilities)
    class(model)<-"mixHMModel"
    model
  }