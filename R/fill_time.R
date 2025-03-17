#' Expand the data.frame with missing time points
#'
#' Base R version of internal `data.table` based `fill_time` function of the 
#' `dynamite` package.
#' 
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#'@noRd 
fill_time <- function(data, time_var, id_var) {
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L, 
    "There must be at least two time points in the data."
  )
  time_ivals <- diff(time)
  time_resolution <- min(time_ivals)
  stopifnot_(
    all(time_ivals %% time_resolution == 0), 
    "Observations must occur at regular time intervals."
  )
  timetable <- table(data[[id_var]], data[[time_var]])
  d <- which(rowSums(timetable > 1) > 0)
  stopifnot_(
    length(d) == 0,
    c("Each time index must correspond to a single observation per group:", 
      x = "{cli::qty(length(d))}ID{?s} {.var {d}} of {.var {id_var}}\n 
    {cli::qty(length(d))}{?has/have} duplicate time points.")
  )
  full_time <- seq(time[1L], time[length(time)], by = time_resolution)
  
  if (sum(timetable) != prod(dim(timetable)) || length(time) != length(full_time)) {
    all_times <- expand.grid(
      time = full_time,
      group = unique(data[[id_var]])
    )
    names(all_times) <- c(time_var, id_var)
    cols <- names(data)
    col_ids <- which(cols %in% c(id_var, time_var))
    idx <- match(
      paste(all_times[[id_var]], all_times[[time_var]]), 
      paste(data[[id_var]], data[[time_var]])
    )
    data <- cbind(
      all_times, 
      data[, -col_ids, drop = FALSE][idx, , drop = FALSE]
    )[, cols]
    rownames(data) <- NULL
  }
  data
}
