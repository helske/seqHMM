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
  type <- try(match.arg(type, c("viterbi", "posterior")), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either {.val viterbi} or {.val posterior}."
  )
  if (type == "viterbi") {
    if (is.null(hp)) {
      hp <- hidden_paths(x)
    } else {
      if (inherits(x, "mhmm")) {
        id <- x$id_variable
        time <- x$time_variable
      } else {
        id <- "id"
        time <- "time"
      }
      cols <- c(id, time, "state", "cluster")
      stopifnot_(
        inherits(hp, "data.table") && all(cols %in% names(hp)),
        "Argument {.arg hp} must be a {.cls data.table} object from 
        {.fun hidden_paths}."
      )
    }
    if (inherits(x, "mnhmm")) id <- x$id_variable else id <- "id"
    clusters <- hp[, .SD[1], by = id, showProgress = FALSE]$cluster
  } else {
    d <- posterior_cluster_probabilities(x)
    clusters <- d[, .SD[which.max(probability)], by = id, 
                  showProgress = FALSE]$cluster
  }
  clusters
}
#' Extract Posterior Cluster Probabilities
#' 
#' @param x An object of class `mhmm` or `mnhmm`.
#' @return a `data.frame` of posterior cluster probabilities for each sequence and 
#' cluster.
#' @export
posterior_cluster_probabilities <- function(x) {
  # avoid CRAN check warnings due to NSE
  probability <- id <- time <- cluster <- NULL
  stopifnot_(
    inherits(x, "mhmm") || inherits(x, "mnhmm"),
    "Argument {.arg x} must be a {.cls mhmm} or {.cls mnhmm} object."
  )
  pp <- posterior_probs(x)
  if (inherits(x, "mhmm")) {
    pp <- pp[time == min(time), list(probability = sum(probability)), 
       by = c("id", "cluster"), showProgress = FALSE]
  } else {
    pp <- pp[time == min(time), list(probability = sum(probability)), 
       by = list(id, cluster), 
       env = list(id = x$id_variable, time = x$time_variable), 
       showProgress = FALSE]
  }
  pp[, cluster := factor(cluster, levels = x$cluster_names)]
  pp[]
}