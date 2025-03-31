#' Transform gammas based on centered QR to original scale
#' @noRd
gamma_std_to_gamma <- function(gamma, R_inv, coef_names, X_mean) {
  non_icpt <- which(coef_names != "(Intercept)")
  icpt <- which(coef_names == "(Intercept)")
  # multichannel or bootstrap samples
  if (is.list(gamma)) {
    gamma <- lapply(
      gamma, gamma_std_to_gamma, 
      R_inv = R_inv, coef_names = coef_names, X_mean = X_mean
    )
  } else {
    if (length(non_icpt) > 0) {
      # pi
      if (is.na(dim(gamma)[3])) {
        for (s in seq_len(nrow(gamma))) {
          gamma[s, non_icpt] <- R_inv %*% gamma[s, non_icpt]
          if (length(icpt) > 0)
            gamma[s, 1] <- gamma[s, 1] - gamma[s, non_icpt] %*% X_mean
        }
      } else {
        # A and B
        for (r in seq_len(dim(gamma)[3])) {
          for (s in seq_len(nrow(gamma))) {
            gamma[s, non_icpt, r] <- R_inv %*% gamma[s, non_icpt, r]
            if (length(icpt) > 0)
              gamma[s, 1, r] <- gamma[s, 1, r] - gamma[s, non_icpt, r] %*% X_mean
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
  non_icpt <- which(coef_names != "(Intercept)")
  icpt <- which(coef_names == "(Intercept)")
  # multichannel or bootstrap samples
  if (is.list(gamma)) {
    gamma <- lapply(
      gamma, gamma_to_gamma_std, 
      R = R, coef_names = coef_names, X_mean = X_mean
    )
  } else {
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
