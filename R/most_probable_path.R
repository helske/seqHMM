#' Most Probable Path of Hidden States of Hidden Markov Model given the
#' Sequence.
#' 
#' Function \code{most_probable_path} computes the most probable path of the
#' hidden states of the hidden Markov model given the observed sequence.
#' 
#' @export
#' @param model Hidden Markov model of class \code{hmm} or
#'  mixture HMM of class \code{mhmm}.
#' 
#' @return List which contains the most probable paths of hidden states (\code{mpp}) 
#'   given the observations and the log-probability (\code{log_prob}). 
#'   For mixture hidden Markov models also the most probable clusters (\code{clusters}) 
#'   and a matrix of cluster probabilities (in columns) by the most probable cluster 
#'   (in rows) (\code{classification_probabilities}). In a case of multiple 
#'   subjects, the most probable path is computed independently for each subject.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' # Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' 
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Starting values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Starting values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' init <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = init
#'   )
#'   
#' # Fitting hidden Markov model 
#' HMM <- fit_hmm(bHMM)
#'   
#' # Computing the most probable paths 
#' mpp <- most_probable_path(HMM$model)$mpp
#'   
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models, \code{\link{build_mhmm}} and 
#'   \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov models, 
#'   and \code{\link{seqIplot}}, \code{\link{ssplot}}, or \code{\link{mssplot}}
#'   for plotting the most probable paths.
#'   

most_probable_path<-function(model){
  
  ll <- logLik(model, partials = TRUE)
  if(inherits(model,"mhmm")){
    fw <- forwardProbs(model)[,model$length_of_sequences,]
    model <- combine_models(model)
    mix<-TRUE
  } else mix <- FALSE
  
  
  if(model$number_of_channels == 1){
    model$observations <- list(model$observations)
    model$emission_matrix <- list(model$emission_matrix)
  }
  
  
  model$initial_probs <- log(model$initial_probs)
  model$transition_matrix <- log(model$transition_matrix)
  
  obsArray<-array(0,c(model$number_of_sequences,model$length_of_sequences,model$number_of_channels))
  for(i in 1:model$number_of_channels){
    obsArray[,,i]<-data.matrix(model$observations[[i]])-1
    obsArray[,,i][obsArray[,,i]>model$number_of_symbols[i]]<-model$number_of_symbols[i]
  }       
  storage.mode(obsArray)<-"integer"
  
  emissionArray<-array(0,c(model$number_of_states,max(model$number_of_symbols)+1,model$number_of_channels))
  for(i in 1:model$number_of_channels)
    emissionArray[,1:model$number_of_symbols[i],i] <- log(model$emission_matrix[[i]])
  
  if(mix){
    out<-viterbix(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray, model$beta, 
      model$X, model$number_of_states_in_clusters)
  } else{
    out<-viterbi(model$transition_matrix, emissionArray, 
      model$initial_probs, obsArray)
  }
  
  
  if(model$number_of_sequences==1){
    mpp<-t(model$state_names[out$q+1])
  }else{
    mpp<-apply(out$q+1,2,function(x) model$state_names[x])
  }
  mpp<-seqdef(mpp,alphabet=model$state_names,
    id=rownames(model$obs[[1]]),
    start=attr(model$obs[[1]],"start"),
    xtstep=attr(model$obs[[1]],"xtstep"))
  
  if(mix==TRUE){
    gr <- sub("^.*?_","",mpp[,1])
    gr <- factor(gr, levels=1:model$number_of_clusters, labels=model$cluster_names)
    
    clP <- vector("list", model$number_of_clusters)
    p <- 0

    for(i in 1:model$number_of_clusters){
      clP[[i]] <- colSums(exp(fw[(p+1):(p+model$number_of_states_in_clusters[i]),] - 
                                rep(ll, each = model$number_of_states_in_clusters[i])))
      p <- p + model$number_of_states_in_clusters[i]
    }
    clProbs <- matrix(NA, nrow = model$number_of_clusters, ncol = model$number_of_clusters)
    rownames(clProbs) <- colnames(clProbs) <- model$cluster_names
    for(i in 1:model$number_of_clusters){
      for(j in 1:model$number_of_clusters){
        clProbs[i,j] <- mean(clP[[j]][gr == model$cluster_names[i]])
      }
    }
    
    list(mpp=mpp, cluster=gr, classification_probabilities = clProbs, log_prob=out$logp)
  }else{
    list(mpp=mpp, log_prob=out$logp)
  }
  
}
