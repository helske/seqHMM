#' @export
#' @rdname print
#' @method print summary.mhmm

print.summary.mhmm <- function(x, digits = 3, ...){
  if(exists("transition_probs", x)){
    cat("Initial probabilities :\n")
    print.listof(x$initial_probs, digits = digits, ...)
    cat("Transition probabilities :\n")
    print.listof(x$transition_probs, digits = digits, ...)
    cat("Emission probabilities :\n")
    if(!is.list(x$emission_probs[[1]])){
      print.listof(x$emission_probs, digits = digits, ...)
    }else{
      for(i in 1:length(x$emission_probs)){
        cat(names(x$emission_probs)[i], ":\n\n")
        print.listof(x$emission_probs[[i]],  digits = digits, ...)
      }
    }
  }
  
  coef_se <- matrix(sqrt(diag(x$vcov)), nrow(x$coefficients))
  coefs <- replicate((ncol(x$coefficients) - 1), 
                     matrix(NA, nrow = nrow(x$coefficients), ncol = 2), simplify = FALSE)
  for(i in 1:length(coefs)){
    coefs[[i]][, 1] <- x$coefficients[, i + 1]
    coefs[[i]][, 2] <- coef_se[, i]
    rownames(coefs[[i]]) <- rownames(x$coefficients)
    colnames(coefs[[i]]) <- c("Estimate", "Std. error")
  }
  cluster_names <- colnames(x$coefficients)
  names(coefs) <- cluster_names[-1]
  cat("Covariate effects :\n")
  cat(cluster_names[1], "is the reference.\n\n")
  print.listof(coefs, print.gap = 2, digits = digits, quote = FALSE, ...)
  
  cat("Log-likelihood:", x$logLik, "  BIC:", x$BIC, "\n\n")
  
  cat("Means of prior cluster probabilities :\n")
  print(colMeans(x$prior_cluster_probabilities), digits = digits, ...)
  cat("\n")
  
  tbl <- table(x$most_probable_cluster)
  cl <- matrix(c(as.character(tbl), 
                 as.character(round(prop.table(tbl), digits = digits))), 
               nrow = 2, byrow = TRUE)
  colnames(cl) <- cluster_names
  rownames(cl) <- c("count", "proportion")
  cat("Most probable clusters :\n")
  print.default(cl, quote = FALSE, print.gap = 2, right = TRUE)
  cat("\n")
  
  cat("Classification table :\n")
  cat("Mean cluster probabilities (in columns) by the most probable cluster (rows)\n\n")
  print(x$classification_table, digits = digits, ...)

}
