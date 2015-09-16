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
#' # Loading a hidden Markov model of the biofam data (hmm object)
#' data(hmm_biofam)
#' 
#' sc <- mc_to_sc(hmm_biofam)
#' 
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models; and \code{\link{hmm_biofam}} for information on 
#'   the model used in the example;.

mc_to_sc<-function(model, combine_missing=TRUE, all_combinations=FALSE){
  
  if(model$n_channels==1)
    return(model)
  
  if(inherits(model, "hmm")){
    
    B<-matrix(0,model$n_states,prod(model$n_symbols))
    
    colnames(B)<-apply(
      expand.grid(lapply(model$emission_matrix,colnames)),                
      1,paste0,collapse="/")
    rownames(B)<-rownames(model$emission_matrix[[1]])
    for(i in 1:model$n_states){
      B[i,]<-apply(expand.grid(lapply(model$emission_matrix,function(x) x[i,])),1,prod)   
    }
    B <- B[, order(colnames(B))]
    
    modelx<-model
    modelx$emission_matrix <- B
    modelx$n_symbols <- ncol(B)
    modelx$n_channels <- as.integer(1)
    modelx$symbol_names <- snames <- colnames(B)
    modelx$channel_names <- "Observations"
    
    modelx$observations<-model$observations[[1]]
    for(i in 2:model$n_channels)
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
    
    cpal <- colorpalette[[modelx$n_symbols]]
    
    if(all_combinations==TRUE){
      modelx$observations<-suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    }else{
      modelx$observations<-suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      modelx$emission_matrix <- modelx$emission_matrix[,colnames(modelx$emission_matrix) %in% alphabet(modelx$observations)==TRUE]
      modelx$symbol_names <- colnames(modelx$emission_matrix)
      modelx$n_symbols <- ncol(modelx$emission_matrix)
    }
    
  # mhmm
  }else{
    
    modelx<-model
    
    B <- vector("list", model$n_clusters)
    for(m in 1:model$n_clusters){
      B[[m]] <- matrix(0,model$n_states[m],prod(model$n_symbols))
      
      colnames(B[[m]])<-apply(
        expand.grid(lapply(model$emission_matrix[[m]],colnames)),                
        1,paste0,collapse="/")
      rownames(B[[m]])<-rownames(model$emission_matrix[[m]][[1]])
      for(i in 1:model$n_states[[m]]){
        B[[m]][i,]<-apply(expand.grid(lapply(model$emission_matrix[[m]],function(x) x[i,])),1,prod)   
      }
      B[[m]] <- B[[m]][, order(colnames(B[[m]]))]

      modelx$emission_matrix[[m]] <- B[[m]]
    }
    
    modelx$n_symbols <- ncol(B[[1]])
    modelx$n_channels <- as.integer(1)
    modelx$symbol_names <- snames <- colnames(B[[1]])
    
    modelx$channel_names <- "Observations"
    
    modelx$observations <- model$observations[[1]]
    for(i in 2:model$n_channels)
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
    
    cpal <- colorpalette[[modelx$n_symbols]]
    
    if(all_combinations==TRUE){
      modelx$observations <- suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    }else{
      modelx$observations<-suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      modelx$emission_matrix <- modelx$emission_matrix[[1]][,colnames(modelx$emission_matrix[[1]]) %in% alphabet(modelx$observations)==TRUE]
      modelx$symbol_names <- colnames(modelx$emission_matrix)
      modelx$n_symbols <- ncol(modelx$emission_matrix)
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

