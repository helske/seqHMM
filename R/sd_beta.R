sd_beta <- function(model){
  sqrt(diag(varcoef(model$beta, model$X, model$number_of_states)))
}