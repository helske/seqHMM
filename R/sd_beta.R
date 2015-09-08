sd_beta <- function(model){
  matrix(c(rep(0,model$number_of_covariates),
    sqrt(diag(varcoef(model$beta, model$X, model$number_of_states)))),
    nrow = model$number_of_covariates, ncol = model$number_of_clusters)
}