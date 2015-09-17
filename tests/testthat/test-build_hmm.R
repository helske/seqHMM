context("Building HMMs")

set.seed(1)

fair <- sample(0:1, size = 50*100, replace = TRUE)
obs <- matrix(0, 50, 100)
obs[fair == 1] <- sample(c("heads", "tails"), size = sum(fair), replace = TRUE)
obs[fair == 0] <- 
  sample(c("heads", "tails"), prob = c(0.7, 0.3), size = sum(fair == 0), replace = TRUE)

emission_matrix <- matrix(c(0.5,0.2,0.5,0.8),2, 2)
transition_matrix <- matrix(c(0.8,0.4,0.2,0.6), 2, 2)
initial_probs <- c(0.5,0.5)

emission_matrix <- matrix(c(0.5,0.2,0.5,0.8),2, 2)
transition_matrix <- matrix(c(0.8,0.3,0.2,0.7), 2, 2)
initial_probs <- c(1,0)

set.seed(1)
sim <- simulate_hmm(n=1,initial_probs, transition_matrix, emission_matrix, 500)


bhmm <- build_hmm(sim$obs, transition_matrix = transition_matrix, 
  emission_matrix = emission_matrix, initial_probs = initial_probs)

hmm <- fit_hmm(bhmm, global = FALSE,local=FALSE)
hmm2 <- fit_hmm(bhmm,hmm$model, em = FALSE, global = FALSE,local=TRUE)

hmm$model
forward_probs(hmm$model)
forward_probs(hmm2$model)

hidden_paths(hmm$model)
hmm2$model
plot(hmm$model)


init_hmm <- initHMM(States = 1:2, Symbols = 1:2, transProbs = hmm2$model$transition_matrix, 
  emissionProbs = hmm2$model$emission_matrix, start = hmm2$model$initial_probs)

bw <- baumWelch(init_hmm, observation = unlist(sim$obs))

# HMM package

init_hmm <- initHMM(States = 1:2, Symbols = 1:2, transProbs = transition_matrix, 
  emissionProbs = emission_matrix, start = initial_probs)

bw <- baumWelch(init_hmm, observation = unlist(sim$obs))


mpp <- hidden_paths(fit$model)
mssplot(fit$model, mpp, plots = "mpp")
test_that("bogus arguments give errors",{  
  expect_error(build_hmm())
})
