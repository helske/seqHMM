#' Convert Coefficients from Stan Output to C++ Input format
#' 
#' @noRd
stan_to_cpp_initial <- function(x, D = 1) {
  if (D > 1) {
    # x is D x S - 1 x K
    K <- dim(x)[3]
    S <- ncol(x) + 1
    z <- matrix(aperm(x, c(2, 1, 3)), D * (S - 1), K)
  } else {
    # x is S - 1 x K
    z <- x
  }
  # output is S - 1 x K or D * (S - 1) x K
  z
}
#' Convert Coefficients from Stan Output to C++ Input format
#' 
#' @noRd
stan_to_cpp_transition <- function(x, D = 1) {
  if (D > 1) {
    # x is D x S x S - 1 x K
    S <- ncol(x)
    K <- dim(x)[4]
    z <- array(0, c(D * (S - 1), K, S))
    for (d in 1:D) {
      z[1:(S - 1) + (d - 1) * (S - 1), ,] <- 
        aperm(array(x[d, , , ], c(S, S - 1, K)), c(2, 3, 1))
    }
  } else {
    # x is S x S - 1 x K
    z <- aperm(x, c(2, 3, 1))
  }
  # output is (S - 1) x K x S or D * (S - 1) x K x S 
  z
}
#' Convert Coefficients from Stan Output to C++ Input format
#' 
#' @noRd
stan_to_cpp_emission <- function(x, D = 1, multichannel = FALSE) {
  if (multichannel) {
    if (D > 1) {
      z <- t(x)
    } else {
      z <- x
    }
  } else {
    if (D > 1) {
      # x is D x S x M - 1 x K
      S <- ncol(x)
      M <- dim(x)[3] + 1
      K <- dim(x)[4]
      z <- array(0, c(D * (M - 1), K, S))
      for (d in 1:D) {
        z[1:(M - 1) + (d - 1) * (M - 1), ,] <- 
          aperm(array(x[d, , , ], c(S, M - 1, K)), c(2, 3, 1))
      }
    } else {
      # x is S x M - 1 x K
      z <- aperm(x, c(2, 3, 1))
    }
    # output is (M - 1) x K x S or D * (M - 1) x K x S 
    z
  }
}