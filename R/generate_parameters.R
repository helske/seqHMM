# if (any(!is.finite(out$theta_tilde))) {
#   warning(paste(
#     "Nonfinite values in samples from normal approximation of the parameters.",
#     "Cannot compute confidence intervals, returning samples for diagnostics.")
#   )
#   list(estimates = out$par[vars], 
#        loglik = ifelse(restarts == 1L, out$value, logliks), 
#        return_code = ifelse(restarts == 1L, out$return_code, return_codes),
#        samples = samples)
# } else {
#   samples <- lapply(
#     as_draws_rvars(out$theta_tilde)[vars],
#     draws_of
#   )
#   cis <- lapply(samples, function(x) {
#     apply(x, 2:length(dim(x)), quantile, probs = c(0.025, 0.975))
#   })