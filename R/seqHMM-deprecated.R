#' Deprecated function(s) in the seqHMM package
#' 
#' These functions are provided for compatibility with older version of
#' the seqHMM package.  They will be eventually completely
#' removed.
#' 
#' @rdname seqHMM-deprecated
#' @name seqHMM-deprecated 
#' @docType package
#' @inheritParams fit_model
#' @export
fit_hmm <- function(model, em_step = TRUE, global_step = FALSE, local_step = FALSE, 
  control_em = list(), control_global = list(), control_local = list(), lb, ub, threads = 1, 
  log_space = FALSE, ...){
  
  .Deprecated("fit_model", package = "seqHMM")
  fit_model(model = model, em_step = em_step, global_step = global_step, local_step = local_step,
    control_em = control_em, control_global = control_global, control_local = control_local,
    lb, ub, threads = threads, log_space = log_space, ...)
}
#' @rdname seqHMM-deprecated
#' @name seqHMM-deprecated 
#' @docType package
#' @inheritParams fit_model
#' @export
fit_mhmm <- function(model, em_step = TRUE, global_step = FALSE, local_step = FALSE, 
  control_em = list(), control_global = list(), control_local = list(), lb, ub, threads = 1, 
  log_space = FALSE, ...){
  
  .Deprecated("fit_model", package = "seqHMM")
  fit_model(model = model, em_step = em_step, global_step = global_step, local_step = local_step,
    control_em = control_em, control_global = control_global, control_local = control_local,
    lb, ub, threads = threads, log_space = log_space, ...)
}

#' @rdname seqHMM-deprecated
#' @name seqHMM-deprecated 
#' @docType package
#' @inheritParams trim_model
#' @export
trim_hmm <- function(model, maxit = 0, return_loglik=FALSE, zerotol=1e-8, verbose = TRUE, ...){
  
  .Deprecated("trim_model", package = "seqHMM")
  trim_model(model = model, maxit = maxit, return_loglik = return_loglik, zerotol = zerotol, 
    verbose = verbose, ...)
}