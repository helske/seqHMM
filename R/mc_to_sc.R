#' Transform Multichannel Hidden Markov Model to Single Channel Representation
#' 
#' Transforms data and parameters of multichannel model to single channel model.
#' Observed states (symbols) are combined and parameters multiplied accross channels.
#' 
#' @export
#' @param model Object of class \code{hmm} or \code{mhmm}.
#' @param combine_missing Controls whether combined states of observations are 
#'   coded missing (NA) if some of the channels include missing information. 
#'   Defaults to \code{TRUE}.
#' @param all_combinations Controls whether all possible combinations of
#'   observed states are included in the single channel representation or only
#'   combinations that are found in the data. Defaults to \code{FALSE}, i.e.
#'   only actual observations are included.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
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
#' ## Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Initial values for emission matrices
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
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities
#' initial_probs <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = initial_probs
#'   )
#' 
#' # Fitting hidden Markov model
#' HMM <- fit_hmm(bHMM)
#' 
#' HMM_SC <- mc_to_sc(HMM$model)
#' 
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models.

mc_to_sc<-function(model, combine_missing=TRUE, all_combinations=FALSE){
  
  if(model$number_of_channels==1)
    return(model)
  
  if(inherits(model, "hmm")){
    
    B<-matrix(0,model$number_of_states,prod(model$number_of_symbols))
    
    colnames(B)<-apply(
      expand.grid(lapply(model$emission_matrix,colnames)),                
      1,paste0,collapse="/")
    rownames(B)<-rownames(model$emission_matrix[[1]])
    for(i in 1:model$number_of_states){
      B[i,]<-apply(expand.grid(lapply(model$emission_matrix,function(x) x[i,])),1,prod)   
    }
    B <- B[, order(colnames(B))]
    
    modelx<-model
    modelx$emission_matrix <- B
    modelx$number_of_symbols <- ncol(B)
    modelx$number_of_channels <- as.integer(1)
    modelx$symbol_names <- snames <- colnames(B)
    modelx$channel_names <- "Observations"
    
    modelx$observations<-model$observations[[1]]
    for(i in 2:model$number_of_channels)
      modelx$observations<-as.data.frame(mapply(paste, modelx$observations,
                                                model$observations[[i]],
                                                USE.NAMES=FALSE,SIMPLIFY=FALSE,
                                                MoreArgs=list(sep="/")))
    names(modelx$observations)<-names(model$observations[[1]])   
    if(combine_missing==TRUE){
      modelx$observations[Reduce("|",
                                 lapply(
                                   model$observations, 
                                   function(x) 
                                     x==attr(model$observations[[1]], "nr") |
                                     x==attr(model$observations[[1]], "void") |
                                     is.na(x)))]<-NA
    }
    
    cpal <- seqHMM::colorpalette[[modelx$number_of_symbols]]
    
    if(all_combinations==TRUE){
      modelx$observations<-suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    }else{
      modelx$observations<-suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      modelx$emission_matrix <- modelx$emission_matrix[,colnames(modelx$emission_matrix) %in% alphabet(modelx$observations)==TRUE]
      modelx$symbol_names <- colnames(modelx$emission_matrix)
      modelx$number_of_symbols <- ncol(modelx$emission_matrix)
    }
    
  # mhmm
  }else{
    
    modelx<-model
    
    B <- vector("list", model$number_of_clusters)
    for(m in 1:model$number_of_clusters){
      B[[m]] <- matrix(0,model$number_of_states_in_clusters[m],prod(model$number_of_symbols))
      
      colnames(B[[m]])<-apply(
        expand.grid(lapply(model$emission_matrix[[m]],colnames)),                
        1,paste0,collapse="/")
      rownames(B[[m]])<-rownames(model$emission_matrix[[m]][[1]])
      for(i in 1:model$number_of_states_in_clusters[[m]]){
        B[[m]][i,]<-apply(expand.grid(lapply(model$emission_matrix[[m]],function(x) x[i,])),1,prod)   
      }
      B[[m]] <- B[[m]][, order(colnames(B[[m]]))]

      modelx$emission_matrix[[m]] <- B[[m]]
    }
    
    modelx$number_of_symbols <- ncol(B[[1]])
    modelx$number_of_channels <- as.integer(1)
    modelx$symbol_names <- snames <- colnames(B[[1]])
    
    modelx$channel_names <- "Observations"
    
    modelx$observations <- model$observations[[1]]
    for(i in 2:model$number_of_channels)
      modelx$observations <- as.data.frame(mapply(paste, modelx$observations,
                                                model$observations[[i]],
                                                USE.NAMES=FALSE,SIMPLIFY=FALSE,
                                                MoreArgs=list(sep="/")))
    names(modelx$observations) <- names(model$observations[[1]])   
    if(combine_missing==TRUE){
      modelx$observations[Reduce("|",
                                 lapply(
                                   model$observations, 
                                   function(x) 
                                     x==attr(model$observations[[1]], "nr") |
                                     x==attr(model$observations[[1]], "void") |
                                     is.na(x)))]<-NA
    }
    
    cpal <- seqHMM::colorpalette[[modelx$number_of_symbols]]
    
    if(all_combinations==TRUE){
      modelx$observations <- suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    }else{
      modelx$observations<-suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      modelx$emission_matrix <- modelx$emission_matrix[[1]][,colnames(modelx$emission_matrix[[1]]) %in% alphabet(modelx$observations)==TRUE]
      modelx$symbol_names <- colnames(modelx$emission_matrix)
      modelx$number_of_symbols <- ncol(modelx$emission_matrix)
    }
  }
  
  attr(modelx$observations, "xtstep") <- attr(model$observations[[1]], "xtstep")
  attr(modelx$observations, "missing.color") <- attr(model$observations[[1]], "missing.color")
  attr(modelx$observations, "nr") <- attr(model$observations[[1]], "nr")
  attr(modelx$observations, "void") <- attr(model$observations[[1]], "void")
  attr(modelx$observations, "missing") <- attr(model$observations[[1]], "missing")
  attr(modelx$observations, "start") <- attr(model$observations[[1]], "start")
  attr(modelx$observations, "cpal") <- cpal[snames %in% modelx$symbol_names]
  modelx
}

