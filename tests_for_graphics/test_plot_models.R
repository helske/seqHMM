# Hidden Markov model

data("hmm_biofam")
plot(hmm_biofam)
plot(hmm_biofam, layout = "vertical", withlegend = "right")
plot(hmm_biofam, )
plot(hmm_biofam, layout = "vertical", pie = FALSE, vertex.label = "names",
  vertex.label.pos = pi/5, vertex.label.family = "serif", loops = TRUE)

data("hmm_mvad")
plot(hmm_mvad)
require("igraph")
plot(hmm_mvad, layout = layout_in_circle, edge.label = NA, legend.prop = 0.25,
  vertex.label.pos = c("bottom", "right", "top", "bottom", "right"),
  withlegend = "left", edge.curved = 0.2)
plot(hmm_mvad, layout = layout_in_circle, label.signif = 2,
  label.max.length = 3, withlegend = FALSE)
plot(hmm_mvad, layout = layout_in_circle, label.signif = 3,
  label.scientific = TRUE, label.max.length = 6, withlegend = FALSE)



# Markov model

data("mvad", package = "TraMineR")

mvad_alphabet <-
  c("employment", "FE", "HE", "joblessness", "school", "training")
mvad_labels <- c("employment", "further education", "higher education",
  "joblessness", "school", "training")
mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet,
  states = mvad_scodes, labels = mvad_labels, xtstep = 6)
attr(mvad_seq, "cpal") <- colorpalette[[6]]
mm_mvad <- build_mm(observations = mvad_seq,
  transition_probs = simulate_transition_probs(6),
  initial_probs = rep(1/6,6))
plot(mm_mvad)
plot(mm_mvad, layout = layout_in_circle, edge.label = NA, withlegend = "bottom")
plot(mm_mvad, layout = layout_in_circle, edge.label = 1:30, withlegend = "right")
plot(mm_mvad, layout = layout_in_circle, edge.label = NA, withlegend = "left")
plot(mm_mvad, layout = layout_in_circle, edge.label = NA, withlegend = "top")


# Mixture hidden Markov model

data("mhmm_biofam")
plot(mhmm_biofam, which.plots = 1)
plot(mhmm_biofam, layout = "vertical", withlegend = "right")
plot(mhmm_biofam, ask = TRUE)
plot(mhmm_biofam, interactive = FALSE, which.plots = 2:3,
  layout = layout_as_star, edge.curved = 0, withlegend = "right")
plot(mhmm_biofam, interactive = FALSE, which.plots = 2:3,
  ncol = 3, layout = "vertical", withlegend = "right")

data("mhmm_mvad")
set.seed(123)
plot(mhmm_mvad, interactive = FALSE, layout = layout_nicely,
  withlegend = "right", edge.curved = 0.2, cex.edge.width = 0.5,
  edge.arrow.size = 0.7, vertex.label.pos = "bottom")


# Mixture Markov model

set.seed(123)
mmm_mvad <- build_mmm(observations = mvad_seq,
  transition_probs = simulate_transition_probs(n_states = 6, n_clusters = 2),
  initial_probs = replicate(2, rep(1/6, 6), simplify = FALSE),
  formula = ~male, data = mvad)
mmm_mvad <- fit_model(mmm_mvad)$model
plot(mmm_mvad)
plot(mmm_mvad, interactive = FALSE, withlegend = "right", layout = layout_as_star)


# Latent class model

set.seed(123)
obs <- seqdef(rbind(
  matrix(sample(letters[1:3], 5000, TRUE, prob = c(0.1, 0.6, 0.3)), 500, 10),
  matrix(sample(letters[1:3], 2000, TRUE, prob = c(0.4, 0.4, 0.2)), 200, 10)))
lcm <- build_lcm(obs, simulate_emission_probs(2, 3))
plot(lcm)
plot(lcm, interactive = FALSE, withlegend = "right")
