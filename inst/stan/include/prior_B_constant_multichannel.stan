// priors for (very weak) regularisation
for (c in 1:C) {
  for (s in 1:S) {
    prior_ += dirichlet_lpdf(B[c, s] | rep_vector(2, M));
  }
}
