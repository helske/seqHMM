#' @export
#' @rdname BIC
BIC.mhmm  <-  function(object, trim = FALSE, ...){

  
  if(trim){
    trimmed <- trim_hmm(object, return_loglik=TRUE, ...)
    object <- trimmed$model
    loglik <- trimmed$loglik
  } else loglik <- logLik(object)
  
  
  -2*loglik + log(object$number_of_sequences*object$length_of_sequences)*
    (sum(unlist(object$initial_probs)>0)+sum(unlist(object$transition_matrix)>0)+
       sum(unlist(object$emission_matrix)>0))
}