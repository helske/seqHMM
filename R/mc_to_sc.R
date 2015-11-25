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
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_model}} for building and
#'   fitting Hidden Markov models; and \code{\link{hmm_biofam}} for information on
#'   the model used in the example;.

mc_to_sc<-function(model, combine_missing=TRUE, all_combinations=FALSE){

  if (model$n_channels == 1){
    return(model)
  }

  if (inherits(model, "hmm")) {

    B <- matrix(0,model$n_states,prod(model$n_symbols))

    colnames(B) <- apply(
      expand.grid(lapply(model$emission_probs,colnames)),
      1,paste0,collapse="/")
    rownames(B) <- rownames(model$emission_probs[[1]])
    for (i in 1:model$n_states) {
      B[i,] <- apply(expand.grid(lapply(model$emission_probs,function(x) x[i,])),1,prod)
    }
    B <- B[, order(colnames(B)), drop = FALSE]

    modelx <- model
    modelx$emission_probs <- B
    modelx$n_symbols <- ncol(B)
    modelx$n_channels <- as.integer(1)
    modelx$symbol_names <- snames <- colnames(B)
    modelx$channel_names <- "Observations"

    modelx$observations<-model$observations[[1]]
    for (i in 2:model$n_channels) {
      modelx$observations<-as.data.frame(mapply(paste, modelx$observations,
                                                model$observations[[i]],
                                                USE.NAMES=FALSE,SIMPLIFY=FALSE,
                                                MoreArgs=list(sep="/")))
    }
    names(modelx$observations) <- names(model$observations[[1]])
    if (combine_missing == TRUE) {
      modelx$observations[Reduce("|",
                                 lapply(
                                   model$observations,
                                   function(x)
                                     x==attr(model$observations[[1]], "nr") |
                                     x==attr(model$observations[[1]], "void") |
                                     is.na(x)))]<-NA
    }

    if (modelx$n_symbols <= 200) {
      cpal <- seqHMM::colorpalette[[modelx$n_symbols]]
    } else {
      cp <- NULL
      k <- 200
      p <- 0
      while(modelx$n_symbols - p > 0){
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      cpal <- cp[1:modelx$n_symbols]
    }


    if (all_combinations == TRUE) {
      modelx$observations <- suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    } else {
      modelx$observations <- suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      modelx$emission_probs <- modelx$emission_probs[, colnames(modelx$emission_probs) %in%
                                                       alphabet(modelx$observations) == TRUE]
      modelx$symbol_names <- colnames(modelx$emission_probs)
      modelx$n_symbols <- ncol(modelx$emission_probs)
    }

    # mhmm
  } else {

    modelx <- model

    B <- vector("list", model$n_clusters)
    for (m in 1:model$n_clusters) {
      B[[m]] <- matrix(0, model$n_states[m], prod(model$n_symbols))

      colnames(B[[m]])<-apply(
        expand.grid(lapply(model$emission_probs[[m]],colnames)),
        1,paste0,collapse="/")
      rownames(B[[m]])<-rownames(model$emission_probs[[m]][[1]])
      for(i in 1:model$n_states[[m]]){
        B[[m]][i,]<-apply(expand.grid(lapply(model$emission_probs[[m]],function(x) x[i,])),1,prod)
      }
      B[[m]] <- B[[m]][, order(colnames(B[[m]])), drop = FALSE]

      modelx$emission_probs[[m]] <- B[[m]]
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
    if (modelx$n_symbols <= 200) {
      cpal <- seqHMM::colorpalette[[modelx$n_symbols]]
    } else {
      cp <- NULL
      k <- 200
      p <- 0
      while(modelx$n_symbols - p > 0){
        cp <- c(cp, seqHMM::colorpalette[[k]])
        p <- p + k
        k <- k - 1
      }
      cpal <- cp[1:modelx$n_symbols]
    }


    if (all_combinations == TRUE) {
      modelx$observations <- suppressWarnings(suppressMessages(seqdef(modelx$observations, alphabet=modelx$symbol_names)))
    } else {
      modelx$observations <- suppressWarnings(suppressMessages((seqdef(modelx$observations))))
      for (m in 1:model$n_clusters){
        modelx$emission_probs[[m]] <- modelx$emission_probs[[m]][, colnames(modelx$emission_probs[[m]]) %in%
                                                         alphabet(modelx$observations) == TRUE]
      }
      modelx$symbol_names <- colnames(modelx$emission_probs[[1]])
      modelx$n_symbols <- ncol(modelx$emission_probs[[1]])
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

