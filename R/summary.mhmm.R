#' Summary method for mixture hidden Markov models
#' 
#' Function \code{summary.mhmm} gives a summary of a mixture hidden Markov model.
#'   
#' 
#' @export
#' @method summary mhmm 
#' @param object Mixture hidden Markov model of class \code{mhmm}.
#' @param parameters Whether or not to print parameters of transition, emission, and 
#' initial probabilities.
#' @param digits Number of decimal places to be used when printing parameters. The 
#' default is 3.
#' @param ... Ignored.
#' 
#' @return \describe{
#'    \item{logLik}{Log-likelihood}
#'    \item{BIC}{Bayesian information criterion}
#'    \item{beta}{Coefficients of covariate parameters}
#'    \item{beta_se}{Standard errors for coefficients}
#'    \item{most_probable_cluster}{The most probable cluster according to posterior probabilities}
#'    \item{prior_cluster_probabilities}{Cluster probabilities (mixing proportions) given by the covariates}
#'    \item{posterior_cluster_probabilities}{Posterior class membership probabilities}
#'    \item{classification_table}{Cluster probabilities (columns) by the most probable cluster (rows)}
#'   }
#'  
#'  
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 5) > 0) | 
#'             (rowSums(bf == 7) > 0 & rowSums(bf == 6) > 0),]
#' children[rownames(bf) %in% rownames(div) & bf == 7] <- "Children"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' married[bf == 7] <- "Divorced"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 & rowSums(bf == 3) == 0 &  
#'           rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0) | 
#'          (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 & rowSums(bf == 3) == 0 &  
#'          rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0),]
#' left[rownames(bf) %in% rownames(wp) & bf == 7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start = 15)
#' marr.seq <- seqdef(married, start = 15)
#' left.seq <- seqdef(left, start = 15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' B2_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.01, 0.99), nrow = 4, ncol = 2, byrow = TRUE)
#'                      
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.29, 0.7, 0.01),
#'                    nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE) 
#' 
#' # Cluster 3
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.01, 0.99), nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow = 6, ncol = 3, byrow = TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.50, 0.50,
#'                     0.01, 0.99,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 6, ncol = 2, byrow = TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                  0,    0.9, 0.07, 0.03, 
#'                  0,      0,  0.9,  0.1, 
#'                  0,      0,    0,    1), 
#'              nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                  0,  0.7,  0.1,   0.1, 0.05, 0.05,
#'                  0,    0, 0.85,  0.01,  0.1, 0.04,
#'                  0,    0,    0,   0.9, 0.05, 0.05,
#'                  0,    0,    0,     0,  0.9,  0.1,
#'                  0,    0,    0,     0,    0,    1), 
#'              nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Birth cohort
#' biofam$cohort <- cut(biofam$birthyr, c(1908, 1935, 1945, 1957))
#' biofam$cohort <- factor(
#'   biofam$cohort, labels=c("1909-1935", "1936-1945", "1946-1957")
#' )
#' 
#' # Build mixture HMM
#' bMHMM <- build_mhmm(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = list(A1,A1,A2),
#'   emission_matrix = list(list(B1_child, B1_marr, B1_left),
#'                         list(B2_child, B2_marr, B2_left), 
#'                         list(B3_child, B3_marr, B3_left)),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   formula = ~ sex + cohort, data = biofam,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Left home")
#'   )
#'   
#' # Fitting hidden Markov model 
#' MHMM <- fit_hmm(bMHMM)
#'   
#' # Model summary
#' summary(MHMM$model)
#'   
#' @seealso \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov models.
#'   

summary.mhmm <- function(object, parameters = FALSE, digits = 3, ...){
  
  ll <- logLik(object, partials = TRUE)
  sum_logLik <- sum(ll)
  BIC <- -2*sum_logLik + log(object$n_sequences*object$length_of_sequences)*
    (sum(unlist(object$initial_probs)>0)+sum(unlist(object$transition_matrix)>0)+
        sum(unlist(object$emission_matrix)>0))
  
  beta_se <- sd_beta(object)
  rownames(beta_se) <- rownames(object$beta)
  colnames(beta_se) <- colnames(object$beta)
  
  fw <- forward_probs(object)[,object$length_of_sequences,]
  
  pr <- exp(object$X%*%object$beta)
  prior_cluster_probabilities <- pr/rowSums(pr)


  posterior_cluster_probabilities <- prior_cluster_probabilities
  p <- 0
  for(i in 1:object$n_clusters){
    posterior_cluster_probabilities[,i] <- colSums(exp(fw[(p+1):(p+object$n_states[i]), , drop = FALSE] - 
                              rep(ll, each = object$n_states[i])))
    p <- p + object$n_states[i]
  }
  most_probable_cluster <- factor(apply(posterior_cluster_probabilities, 1, which.max), 
    levels = 1:object$n_clusters, labels = object$cluster_names)
  

  clProbs <- matrix(NA, nrow = object$n_clusters, ncol = object$n_clusters)
  rownames(clProbs) <- colnames(clProbs) <- object$cluster_names
  for(i in 1:object$n_clusters){
    for(j in 1:object$n_clusters){
      clProbs[i,j] <- mean(posterior_cluster_probabilities[most_probable_cluster == object$cluster_names[i], j])
    }
  }
  
  if(!parameters){
    summary_mhmm <- list(
      logLik = sum_logLik, BIC = BIC, most_probable_cluster = most_probable_cluster, 
      beta = object$beta, beta_se = beta_se,
      prior_cluster_probabilities = prior_cluster_probabilities, 
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs,
      digits = digits
    )
  }else{
    summary_mhmm <- list(
      transition_matrix = object$transition_matrix,
      emission_matrix = object$emission_matrix,
      initial_probs = object$initial_probs,
      logLik = sum_logLik, BIC = BIC, most_probable_cluster = most_probable_cluster, 
      beta = object$beta, beta_se = beta_se,
      prior_cluster_probabilities = prior_cluster_probabilities, 
      posterior_cluster_probabilities = posterior_cluster_probabilities,
      classification_table = clProbs,
      digits = digits,
      model = object
    )
  }
  class(summary_mhmm) <- "summary.mhmm"
  summary_mhmm
}