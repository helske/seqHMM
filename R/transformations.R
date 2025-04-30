#' (Regularized) Inverse of softmax(Q*eta)
#' 
#' @noRd
p_to_eta <-function(x) {
  Q <- create_Q(length(x))
  x <- pmin(pmax(x, 1e-6), 1-1e-6)
  x <- x / sum(x)
  log_x <- log(x)
  t(Q) %*% (log_x - mean(log_x))
}

#' Transform gammas based on centered QR to original scale
#' @noRd
gamma_std_to_gamma <- function(gamma, R_inv, coef_names, X_mean) {
  
  # multichannel, mixture or bootstrap samples
  if (is.list(gamma)) {
    # multichannel
    if (!is.list(gamma[[1]]) && is.list(coef_names)) {
      gamma <- lapply(
        seq_along(gamma), 
        \(i) gamma_std_to_gamma(
          gamma[[i]], R_inv[[i]], coef_names[[i]], X_mean[[i]]
        )
      )
      return(gamma)
    }
    # mixture or bootstrap
    gamma <- lapply(
      seq_along(gamma), 
      \(i) gamma_std_to_gamma(gamma[[i]], R_inv, coef_names, X_mean)
    )
  } else {
    non_icpt <- which(coef_names != "(Intercept)")
    icpt <- which(coef_names == "(Intercept)")
    if (length(non_icpt) > 0) {
      # pi
      if (is.na(dim(gamma)[3])) {
        for (s in seq_len(nrow(gamma))) {
          gamma[s, non_icpt] <- R_inv %*% gamma[s, non_icpt]
          if (length(icpt) > 0) {
            gamma[s, 1] <- gamma[s, 1] - gamma[s, non_icpt] %*% X_mean
          }
        }
      } else {
        # A and B
        for (r in seq_len(dim(gamma)[3])) {
          for (s in seq_len(nrow(gamma))) {
            gamma[s, non_icpt, r] <- R_inv %*% gamma[s, non_icpt, r]
            if (length(icpt) > 0) {
              gamma[s, 1, r] <- gamma[s, 1, r] - gamma[s, non_icpt, r] %*% X_mean
            }
          }
        }
      }
    }
  }
  gamma
}
#' Transform natural gammas to centered QR scale
#' @noRd
gamma_to_gamma_std <- function(gamma, R, coef_names, X_mean) {
  # multichannel, mixture or bootstrap samples
  if (is.list(gamma)) {
    # multichannel
    if (!is.list(gamma[[1]]) && is.list(coef_names)) {
      gamma <- lapply(
        seq_along(gamma), 
        \(i) gamma_to_gamma_std(
          gamma[[i]], R[[i]], coef_names[[i]], X_mean[[i]]
        )
      )
      return(gamma)
    }
    # mixture or bootstrap
    gamma <- lapply(
      seq_along(gamma), 
      \(i) gamma_to_gamma_std(gamma[[i]], R, coef_names, X_mean)
    )
  } else {
    non_icpt <- which(coef_names != "(Intercept)")
    icpt <- which(coef_names == "(Intercept)")
    if (length(non_icpt) > 0) {
      # pi
      if (is.na(dim(gamma)[3])) {
        for (s in seq_len(nrow(gamma))) {
          if (length(icpt) > 0) {
            gamma[s, 1] <- gamma[s, 1] + gamma[s, non_icpt] %*% X_mean
          }
          gamma[s, non_icpt] <- R %*% gamma[s, non_icpt]
        }
      } else {
        # A and B
        for (r in seq_len(dim(gamma)[3])) {
          for (s in seq_len(nrow(gamma))) {
            if (length(icpt) > 0) {
              gamma[s, 1, r] <- gamma[s, 1, r] + gamma[s, non_icpt, r] %*% X_mean
            }
            gamma[s, non_icpt, r] <- R %*% gamma[s, non_icpt, r]
          }
        }
      }
    }
  }
  gamma
}