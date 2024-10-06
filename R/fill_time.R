#' Expand the data.frame with missing time points
#'
#' Base R version of internal `data.table` based `fill_time` function of the 
#' `dynamite` package.
#' 
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#'@noRd 
fill_time <- function(data, id_var, time_var) {
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L, 
    "There must be at least two time points in the data."
  )
  time_ivals <- diff(time)
  time_resolution <- min(time_ivals)
  full_time <- seq(time[1L], time[length(time)], by = time_resolution)
  ids <- unique(data[[id_var]])
  n_group <- length(ids)
  time_duplicated <- logical(n_group)
  time_missing <- logical(n_group)
  group_bounds <- c(0L, cumsum(rle(data[[id_var]])$lengths))
  for (i in seq_len(n_group)) {
    idx_group <- seq(group_bounds[i] + 1L, group_bounds[i + 1L])
    sub <- data[idx_group, ]
    time_duplicated[i] <- any(duplicated(sub[[time_var]]))
    time_missing[i] <- !identical(sort(sub[[time_var]]), full_time)
  }
  d <- which(time_duplicated)
  stopifnot_(
    all(!time_duplicated), 
    c("Each time index must correspond to a single observation per group:", 
      x = "{cli::qty(length(d))}ID{?s} {.var {d}} of {.var {id_var}}\n 
    {cli::qty(length(d))}{?has/have} duplicate time points."))
  stopifnot_(
    all(time_ivals[!is.na(time_ivals)]%%time_resolution == 0), 
    "Observations must occur at regular time intervals."
  )
  if (any(time_missing)) {
    all_times <- expand.grid(
      group = unique(data[[id_var]]),
      time = full_time
    )
    names(all_times) <- c(id_var, time_var)
    
    data <- merge(all_times, data, by = c(id_var, time_var), all.x = TRUE)
  }
  data
}
