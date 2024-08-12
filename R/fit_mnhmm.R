#' Build a Mixture Hidden Markov Model with Covariates
#'
#' @noRd
estimate_mnhmm <- function(
    observations, number_of_states, number_of_mixtures, 
    initial_formula, 
    transition_formula, emission_formula, mixture_formula, 
    data = NULL, 
    init_data = NULL, restarts = 1L, threads = 1L, ...) {
  
  multichannel <- is_multichannel(observations)
  # Single channel but observations is a list
  if (is.list(observations) && !inherits(observations, "stslist") && length(observations) == 1) {
    observations <- observations[[1]]
    multichannel <- FALSE
  }
  n_channels <- ifelse(multichannel, length(observations), 1L)
  if(n_channels > 1L) stop("Currently only single-channel sequences are supported.")
  
  N <- nrow(observations)
  T <- ncol(observations)
  M <- length(alphabet(observations))
  observations <- t(sapply(observations, as.integer))
  observations[observations > M] <- 0L # missing values to zero
  
  vars <- c("alpha_m", "beta_m", "alpha_i", "beta_i", "beta_s", "beta_o", "rho", "A", "B")
  # this could be optimized so to avoid duplicate variables 
  # (i.e. use one common X in Stan with suitable subsetting of columns)
  if (is.null(init_data)) {
    X_i <- X_m <- matrix(0, N, 0L)
  } else {
    X_i <- model.matrix.lm(initial_formula, data = init_data, na.action = na.pass)
    X_i[is.na(X_i)] <- 0
    if (ncol(X_i) < 1L) {
      stop("Initial state probability formula should contain at least one term.")
    }
    cnames <- colnames(X_i)
    if ("(Intercept)" %in% cnames) {
      X_i <- X_i[, cnames != "(Intercept)", drop = FALSE]
    }
    X_m <- model.matrix.lm(mixture_formula, data = init_data, na.action = na.pass)
    X_m[is.na(X_m)] <- 0
    if (ncol(X_m) < 1L) {
      stop("Initial state probability formula should contain at least one term.")
    }
    cnames <- colnames(X_m)
    if ("(Intercept)" %in% cnames) {
      X_m <- X_m[, cnames != "(Intercept)", drop = FALSE]
    }
  }
  K_i <- ncol(X_i)
  X_i <- t(X_i)
  K_m <- ncol(X_m)
  X_m <- t(X_m)
  if (K_i == 0L) vars <- vars[-4L]
  if (K_m == 0L) vars <- vars[-2L]
  if (is.null(data)) {
    X_s <- X_o <- array(1, c(T, 1, N))
    K_s <- K_o <- 1L
  } else {
    X_s <- model.matrix.lm(transition_formula, data = data, na.action = na.pass)
    X_s[is.na(X_s)] <- 0
    if (ncol(X_s) < 1L) {
      stop("State transition probability formula should contain at least one term.")
    }
    K_s <- ncol(X_s)
    X_s <- t(X_s)
    dim(X_s) <- c(T, K_s, N)
    
    X_o <- model.matrix.lm(transition_formula, data = data, na.action = na.pass)
    X_o[is.na(X_o)] <- 0
    cnames <- colnames(X_o)
    if (ncol(X_o) < 1L) {
      stop("Emission probability formula should at least one term.")
    }
    K_o <- ncol(X_o)
    X_o <- t(X_o)
    dim(X_o) <- c(T, K_o, N)
  }
  
  if (restarts > 1L) {
    if (threads > 1L) {
      plan(multisession, workers = threads)
    } else {
      plan(sequential)
    }
    out <- future_lapply(seq_len(restarts), function(i) {
      optimizing(
        stanmodels$model, init = "random",
        data = list(
          N = N,
          T = T, 
          M = M,
          S = as.integer(number_of_states),
          D = as.integer(number_of_mixtures),
          K_m = K_m,
          K_i = K_i,
          K_s = K_s,
          K_o = K_o,
          X_m = X_m,
          X_i = X_i,
          X_s = X_s,
          X_o = X_o,
          obs = observations),
        ...)[c("par", "value", "return_code")]
    },
    future.seed = TRUE)
    logliks <- unlist(lapply(out, "[[", "value"))
    return_codes <- unlist(lapply(out, "[[", "return_code"))
    successful <- which(return_codes == 0)
    optimum <- successful[which.max(logliks[successful])]
    init <- as.list(out[[optimum]]$par)
  } else {
    init <- "random"
  }
  
  out <- optimizing(
    stanmodels$model, 
    data = list(
      N = N,
      T = T, 
      M = M,
      S = as.integer(number_of_states),
      D = as.integer(number_of_mixtures),
      K_m = K_m,
      K_i = K_i,
      K_s = K_s,
      K_o = K_o,
      X_m = X_m,
      X_i = X_i,
      X_s = X_s,
      X_o = X_o,X_o = X_o,
      obs = observations), 
    as_vector = FALSE,
    draws = 500,
    init = init,
    ...)[c("par", "value", "return_code", "theta_tilde")]
  
  if (any(!is.finite(out$theta_tilde))) {
    warning(paste(
      "Nonfinite values in samples from normal approximation of the parameters.",
      "Cannot compute confidence intervals, returning samples for diagnostics.")
    )
    list(estimates = out$par[vars], 
         loglik = ifelse(restarts == 1L, out$value, logliks), 
         return_code = ifelse(restarts == 1L, out$return_code, return_codes),
         samples = samples)
  } else {
    samples <- lapply(
      as_draws_rvars(out$theta_tilde)[vars],
      draws_of
    )
    cis <- lapply(samples, function(x) {
      apply(x, 2:length(dim(x)), quantile, probs = c(0.025, 0.975))
    })
    
    list(estimates = out$par[vars], cis = cis, 
         loglik = ifelse(restarts == 1L, out$value, logliks), 
         return_code = ifelse(restarts == 1L, out$return_code, return_codes))
  }
}


