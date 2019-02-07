#' Estimate Regression Coefficients of Mixture Hidden
#' Markov Models
#'
#' Function \code{estimate_coef} estimates the regression coefficients of mixture hidden
#' Markov models and its restricted variants while keeping other parameters fixed.
#'
#' @export
#' @param model An object of class \code{hmm} or \code{mhmm}.
#' @param threads Number of threads to use in parallel computing. The default is 1.
estimate_coef <- function(model, threads = 1){
  
  
  if (!inherits(model, c("hmm", "mhmm"))){
    stop("Argument model must be an object of class 'hmm' or 'mhmm.")
  }
  
  if (threads < 1) {
    stop("Argument threads must be a positive integer.")
  }
  
  df <- attr(model, "df")
  nobs <- attr(model, "nobs")
  original_model <- model
  model <- combine_models(model)
  
  if (model$n_channels == 1) {
    model$observations <- list(model$observations)
    model$emission_probs <- list(model$emission_probs)
  }
  
  obsArray <- array(0, c(model$n_sequences, model$length_of_sequences,
    model$n_channels))
  for(i in 1:model$n_channels) {
    obsArray[,,i] <- sapply(model$observations[[i]], as.integer) - 1L 
    obsArray[,,i][obsArray[,,i] > model$n_symbols[i]] <- model$n_symbols[i]
  }
  obsArray <- aperm(obsArray)
  
  emissionArray <- array(1,c(model$n_states,max(model$n_symbols)+1,model$n_channels))
  for(i in 1:model$n_channels)
    emissionArray[,1:model$n_symbols[i],i] <- model$emission_probs[[i]]
  
  em.con <- list(print_level = 0, maxeval = 1000, reltol = 1e-10)
  
  res <- estimate_coefs(model$transition_probs, emissionArray, model$initial_probs, obsArray,
    model$n_symbols, model$coefficients, model$X, model$n_states_in_clusters,
    em.con$maxeval, em.con$reltol,em.con$print_level, threads)
  
  
  if (res$error != 0) {
    err_msg <- switch(res$error,
      "Scaling factors contain non-finite values.",
      "Backward probabilities contain non-finite values.",
      "Initial values of coefficients of covariates gives non-finite cluster probabilities.",
      "Estimation of coefficients of covariates failed due to singular Hessian.",
      "Estimation of coefficients of covariates failed due to non-finite cluster probabilities.",
      "Non-finite log-likelihood")
    stop(paste("EM algorithm failed:", err_msg))
  }
  
  original_model$coefficients[] <- res$coefficients
  original_model
}
