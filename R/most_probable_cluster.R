#' Extract Most Probable Cluster for Each Sequence
#' 
#' @param x An object of class `mhmm` or `mnhmm`.
#' @param type A character string specifying the method to use. Either
#' `"viterbi"` (default) or `"posterior"`. Former uses the most probable hidden 
#' path to determine the cluster membership for each sequence, while the latter 
#' finds the cluster which has the largest sum of posterior probabilities of 
#' states of that cluster. 
#' @param hp An output from [hidden_paths()] function. Only used in case of 
#' `type = "viterbi"`. If missing, hidden paths will be computed using `x`.
#' @return A vector containing the most probable cluster for each sequence.
#' @export
most_probable_cluster <- function(x, type = "viterbi", hp = NULL) {
  # avoid CRAN check warnings due to NSE
  probability <- NULL
  stopifnot_(
    inherits(x, "mhmm") || inherits(x, "mnhmm"),
    "Argument {.arg x} must be a {.cls mhmm} or {.cls mnhmm} object."
  )
  type <- match.arg(type, c("viterbi", "posterior"))
  if (type == "viterbi") {
    if (is.null(hp)) hp <- hidden_paths(x)
    if (inherits(x, "mnhmm")) id <- x$id_variable else id <- "id"
    clusters <- hp[, .SD[1], by = id]$cluster
  } else {
    d <- posterior_cluster_probabilities(x)
    clusters <- d[, .SD[which.max(probability)], by = id]$cluster
  }
  clusters
}
#' Extract Posterior Cluster Probabilities
#' 
#' @param x An object of class `mhmm` or `mnhmm`.
#' @return matrix of posterior cluster probabilities for each sequence and 
#' cluster.
#' @export
posterior_cluster_probabilities <- function(x) {
  # avoid CRAN check warnings due to NSE
  probability <- id <- time <- cluster <- NULL
  stopifnot_(
    mhmm <- inherits(x, "mhmm") || inherits(x, "mnhmm"),
    "Argument {.arg x} must be a {.cls mhmm} or {.cls mnhmm} object."
  )
  pp <- posterior_probs(x)[time == min(time), ]
  if (mhmm) {
    pp[, list(probability = sum(probability)), by = list(id, cluster)]
  } else {
    pp[, list(probability = sum(probability)), by = c(x$id_variable, "cluster")]
  }
}