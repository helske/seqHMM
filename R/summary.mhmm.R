#' Summary method for mixture hidden Markov models
#' 
#' Function \code{summary.mhmm} gives a summary of a mixture hidden Markov model.
#'   
#' 
#' @export
#' @param model Hidden Markov model of class \code{hmm} or
#'  mixture HMM of class \code{mhmm}.
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
#' # Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' 
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Starting values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Starting values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' init <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = init
#'   )
#'   
#' # Fitting hidden Markov model 
#' HMM <- fit_hmm(bHMM)
#'   
#' # Computing the most probable paths 
#' mpp <- hidden_paths(HMM$model)$mpp
#'   
#' @seealso \code{\link{fit_mhmm}} for building and fitting mixture hidden Markov models.
#'   

summary.mhmm <- function(model, parameters = FALSE){
  
  logLik = logLik(model)
  
  ll <- logLik(model, partials = TRUE)
  
  BIC = BIC(model)
  
  beta_se = sd_beta(model)
  
  fw <- forward_probs(model)[,model$length_of_sequences,]
  
  pr <- exp(model$X%*%model$beta)
  cluster_probabilities <- pr/rowSums(pr)

  
  gr <- sub("^.*?_","",mpp[,1])
  gr <- factor(gr, levels=1:model$number_of_clusters, labels=model$cluster_names)
  
  clP <- vector("list", model$number_of_clusters)
  p <- 0
  
  for(i in 1:model$number_of_clusters){
    clP[[i]] <- colSums(exp(fw[(p+1):(p+model$number_of_states[i]), , drop = FALSE] - 
                              rep(ll, each = model$number_of_states[i])))
    p <- p + model$number_of_states[i]
  }
  clProbs <- matrix(NA, nrow = model$number_of_clusters, ncol = model$number_of_clusters)
  rownames(clProbs) <- colnames(clProbs) <- model$cluster_names
  for(i in 1:model$number_of_clusters){
    for(j in 1:model$number_of_clusters){
      clProbs[i,j] <- mean(clP[[j]][gr == model$cluster_names[i]])
    }
  }
  
  list(logLik = logLik, BIC = BIC, most_probable_cluster = gr, 
       beta = model$beta, beta_se = beta_se,
       prior_cluster_probabilities = cluster_probabilities, 
       posterior_cluster_probabilities = clP,
       classification_table = clProbs
       )
}