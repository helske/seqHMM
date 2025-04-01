library(seqHMM)
data("biofam3c")

# Building sequence objects
marr_seq <- seqdef(
  biofam3c$married,
  start = 15,
  alphabet = c("single", "married", "divorced"),
  cpal = c("violetred2", "darkgoldenrod2", "darkmagenta"
  )
)
child_seq <- seqdef(
  biofam3c$children,
  start = 15,
  alphabet = c("childless", "children"),
  cpal = c("darkseagreen1", "coral3"
  )
)
left_seq <- seqdef(
  biofam3c$left,
  start = 15,
  alphabet = c("with parents", "left home"),
  cpal = c("lightblue", "red3"
  )
)

## Starting values for emission probabilities
# Cluster 1
B1_marr <- matrix(
  c(0.8, 0.1, 0.1, # High probability for single
    0.8, 0.1, 0.1,
    0.3, 0.6, 0.1, # High probability for married
    0.3, 0.3, 0.4), # High probability for divorced
  nrow = 4, ncol = 3, byrow = TRUE)

B1_child <- matrix(
  c(0.9, 0.1, # High probability for childless
    0.9, 0.1,
    0.9, 0.1,
    0.9, 0.1),
  nrow = 4, ncol = 2, byrow = TRUE)

B1_left <- matrix(
  c(0.9, 0.1, # High probability for living with parents
    0.1, 0.9, # High probability for having left home
    0.1, 0.9,
    0.1, 0.9),
  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 2

B2_marr <- matrix(
  c(0.8, 0.1, 0.1, # High probability for single
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1, # High probability for married
    0.7, 0.2, 0.1),
  nrow = 4, ncol = 3, byrow = TRUE)

B2_child <- matrix(
  c(0.9, 0.1, # High probability for childless
    0.9, 0.1,
    0.9, 0.1,
    0.1, 0.9),
  nrow = 4, ncol = 2, byrow = TRUE)

B2_left <- matrix(
  c(0.9, 0.1, # High probability for living with parents
    0.1, 0.9,
    0.1, 0.9,
    0.1, 0.9),
  nrow = 4, ncol = 2, byrow = TRUE)

# Cluster 3
B3_marr <- matrix(
  c(0.8, 0.1, 0.1, # High probability for single
    0.8, 0.1, 0.1,
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1, # High probability for married
    0.3, 0.4, 0.3,
    0.1, 0.1, 0.8), # High probability for divorced
  nrow = 6, ncol = 3, byrow = TRUE)

B3_child <- matrix(
  c(0.9, 0.1, # High probability for childless
    0.9, 0.1,
    0.5, 0.5,
    0.5, 0.5,
    0.5, 0.5,
    0.1, 0.9),
  nrow = 6, ncol = 2, byrow = TRUE)


B3_left <- matrix(
  c(0.9, 0.1, # High probability for living with parents
    0.1, 0.9,
    0.5, 0.5,
    0.5, 0.5,
    0.1, 0.9,
    0.1, 0.9),
  nrow = 6, ncol = 2, byrow = TRUE)

# Starting values for transition matrices
A1 <- matrix(
  c(0.80, 0.16, 0.03, 0.01,
    0,    0.90, 0.07, 0.03,
    0,    0,    0.90, 0.10,
    0,    0,    0,       1),
  nrow = 4, ncol = 4, byrow = TRUE)

A2 <- matrix(
  c(0.80, 0.10, 0.05, 0.03, 0.01, 0.01,
    0,    0.70, 0.10, 0.10, 0.05, 0.05,
    0,    0,    0.85, 0.01, 0.10, 0.04,
    0,    0,    0,    0.90, 0.05, 0.05,
    0,    0,    0,    0,    0.90, 0.10,
    0,    0,    0,    0,    0,       1),
  nrow = 6, ncol = 6, byrow = TRUE)

# Starting values for initial state probabilities
initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)

# Birth cohort
biofam3c$covariates$cohort <- factor(
  cut(biofam3c$covariates$birthyr,
      c(1908, 1935, 1945, 1957), 
      c("1909-1935", "1936-1945", "1946-1957")
  )
)

# Build mixture HMM
init_mhmm_bf <- build_mhmm(
  observations = list(marr_seq, child_seq, left_seq),
  initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
  transition_probs = list(A1, A1, A2),
  emission_probs = list(
    list(B1_marr, B1_child, B1_left),
    list(B2_marr, B2_child, B2_left),
    list(B3_marr, B3_child, B3_left)
  ),
  formula = ~sex + cohort, data = biofam3c$covariates,
  channel_names = c("Marriage", "Parenthood", "Residence"))

# Fitting the model
mhmm_biofam <- fit_model(init_mhmm_bf)$model

usethis::use_data(mhmm_biofam)
