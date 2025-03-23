#' Expand the data.frame with missing time points
#'
#' Slightly modified version of the internal `fill_time` function of the 
#' `dynamite` package.
#' 
#' @references Tikka S, Helske J (2024). “dynamite: An R Package for Dynamic 
#' Multivariate Panel Models.” doi:10.48550/arXiv.2302.01607.
#'@noRd 
fill_time <- function(data, id_var, time_var) {
  times <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(times) > 1L, 
    "There must be at least two time points in the data."
  )
  time_ivals <- diff(times)
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
  full_time <- seq(times[1L], times[length(times)], by = time_resolution)
  
  if (sum(timetable) != prod(dim(timetable)) || length(times) != length(full_time)) {
    data_names <- names(data)
    full_data_template <- as.data.table(
      expand.grid(
        group = unique(data[[id_var]]), 
        time = full_time,
        stringsAsFactors = FALSE)
    )
    names(full_data_template) <- c(id_var, time_var)
    data <- merge.data.table(
      full_data_template, data, by = c(id_var, time_var), all.x = TRUE
    )
    setcolorder(data, data_names)
  }
  data
}
