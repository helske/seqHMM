#' Deprecated function(s) in the seqHMM package
#' 
#' These functions still work but will be removed (defunct) in the next version 
#' of `seqHMM`.
#' 
#' * `ssplot`, `ssp`, `mssplot`, `plot.ssp`. Use [stacked_sequence_plot()] instead.
#' * `gridplot` Use [stacked_sequence_plot()] and `patchwork` package instead.
#' 
#' @name seqHMM-deprecated
NULL

#' Defunct function(s) in the seqHMM package
#'
#' These functions are no longer available in the seqHMM package.
#' 
#' @name seqHMM-defunct
fit_hmm <- function(
    model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
    control_em = list(), control_global = list(), control_local = list(), 
    lb, ub, threads = 1, log_space = FALSE, ...) {
  .Defunct("fit_model", package = "seqHMM")
}
#' @name seqHMM-defunct
fit_mhmm <- function(
    model, em_step = TRUE, global_step = FALSE, local_step = FALSE,
    control_em = list(), control_global = list(), control_local = list(), 
    lb, ub, threads = 1, log_space = FALSE, ...) {
  .Defunct("fit_model", package = "seqHMM")
}
#' @name seqHMM-defunct
trim_hmm <- function(
    model, maxit = 0, return_loglik = FALSE, zerotol = 1e-8, 
    verbose = TRUE, ...) {
  .Defunct("trim_model", package = "seqHMM")
}
