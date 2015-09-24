#' @export
#' @rdname print
#' @method print summary.mhmm

print.summary.mhmm <- function(x, digits = 3, ...){
  if(exists("transition_matrix", x)){
    if(x$model$n_channels == 1){
      cat("Initial probabilities :\n\n")
      print.listof(x$initial_probs, digits = digits, ...)
      cat("Transition probabilities :\n\n")
      print.listof(x$transition_matrix, digits = digits, ...)
      cat("Emission probabilities :\n\n")
      print.listof(x$emission_matrix, digits = digits, ...)
      cat("\n")
    }else{
      cat("Initial probabilities :\n\n")
      print.listof(x$initial_probs, digits = digits, ...)
      cat("Transition matrix :\n\n")
      print.listof(x$transition_matrix, digits = digits, ...)
      cat("Emission matrix :\n\n")
      for(i in 1:length(x$emission_matrix)){
        cat(names(x$emission_matrix)[i], ":\n\n")
        print.listof(x$emission_matrix[[i]],  digits = digits, ...)
      }
      cat("\n")
    }
  }
  
  coefs <- replicate((x$model$n_clusters - 1), 
                     matrix(NA, nrow = nrow(x$coefficients), ncol = 2), simplify = FALSE)
  for(i in 1:length(coefs)){
    coefs[[i]][, 1] <- signif(x$coefficients[, i+1], digits = digits, ...)
    coefs[[i]][, 2] <- signif(x$coef_se[, i + 1], digits = digits, ...)
    rownames(coefs[[i]]) <- rownames(x$coefficients)
    colnames(coefs[[i]]) <- c("Estimate", "Std. error")
  }
  names(coefs) <- x$model$cluster_names[-1]
  cat("Covariate effects :\n")
  cat(x$model$cluster_names[1], "is the reference.\n\n")
  print.listof(coefs, print.gap = 2, digits = digits, ...)
  cat("\n")
  
  cat("Log-likelihood:", x$logLik, "  BIC:", x$BIC, "\n\n")
  
  cat("Means of prior cluster probabilities :\n")
  print(colMeans(x$prior_cluster_probabilities), digits = digits, ...)
  cat("\n\n")
  
  tbl <- table(x$most_probable_cluster)
  cl <- matrix(c(as.character(tbl), 
                 as.character(round(prop.table(tbl), digits = digits))), 
               nrow = 2, byrow = TRUE)
  colnames(cl) <- x$model$cluster_names
  rownames(cl) <- c("count", "proportion")
  cat("Most probable clusters :\n")
  print.default(cl, quote = FALSE, print.gap = 2, right = TRUE)
  cat("\n\n")
  
  cat("Classification table :\n")
  cat("Mean cluster probabilities (in columns) by the most probable cluster (rows)\n\n")
  print(x$classification_table, digits = digits, ...)

}