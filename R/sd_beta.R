sd_coefs <- function(model){
  matrix(c(rep(0,model$n_covariates),
    sqrt(diag(varcoef(model$coefficients, model$X, model$n_states)))),
    nrow = model$n_covariates, ncol = model$n_clusters)
}