#' @export
#' @rdname print
#' @method print summary.mhmm

print.summary.mhmm <- function(x, ...){
  if(exists("transition_matrix", x)){
    if(x$model$n_channels == 1){
      cat("Transition matrix :\n\n")
      print.listof(lapply(x$transition_matrix, round, digits = x$digits))
      cat("Emission matrix :\n\n")
      print.listof(lapply(x$emission_matrix, round, digits = x$digits))
      cat("Initial probabilities :\n\n")
      print.listof(lapply(x$initial_probs, round, digits = x$digits))
      cat("\n\n\n")
    }else{
      cat("Transition matrix :\n\n")
      print.listof(lapply(x$transition_matrix, round, digits = x$digits))
      cat("Emission matrix :\n\n")
      for(i in 1:length(x$emission_matrix)){
        cat(names(x$emission_matrix)[i], ":\n\n")
        print.listof(lapply(x$emission_matrix[[i]], round, digits = x$digits))
      }
      cat("Initial probabilities :\n\n")
      print.listof(lapply(x$initial_probs, round, digits = x$digits))
      cat("\n\n\n")
    }
  }
  
  coeff <- replicate((x$model$n_clusters - 1), 
                     matrix(NA, nrow = nrow(x$beta), ncol = 2), simplify = FALSE)
  for(i in 1:length(coeff)){
    coeff[[i]][, 1] <- signif(x$beta[, i+1], digits = x$digits)
    coeff[[i]][, 2] <- signif(x$beta_se[, i+1], digits = x$digits)
    rownames(coeff[[i]]) <- rownames(x$beta)
    colnames(coeff[[i]]) <- c("Estimate", "Std. error")
  }
  names(coeff) <- x$model$cluster_names[-1]
  cat("Covariate effects :\n\n")
  cat(x$model$cluster_names[1], "is the reference.\n\n")
  print.listof(coeff, print.gap = 2)
  cat("\n")
  
  cat("Log-likelihood :\n\n")
  cat(x$logLik)
  cat("\n\n\n")
  
  cat("BIC :\n\n")
  cat(x$BIC)
  cat("\n\n\n")
  
  cat("Means of prior cluster probabilities :\n\n")
  print(round(colMeans(x$prior_cluster_probabilities), digits = x$digits))
  cat("\n\n")
  
  tbl <- table(x$most_probable_cluster)
  cl <- matrix(c(as.character(tbl), 
                 as.character(round(prop.table(tbl)*100, digits = x$digits))), 
               nrow = 2, byrow = TRUE)
  colnames(cl) <- x$model$cluster_names
  rownames(cl) <- c("n", "%")
  cat("Most probable clusters :\n\n")
  print.default(cl, quote = FALSE, print.gap = 2, right = TRUE)
  cat("\n\n")
  
  cat("Classification table :\n")
  cat("Mean cluster probabilities (in columns) by the most probable cluster (rows)\n\n")
  print(round(x$classification_table, digits = x$digits))

}