// priors for (very weak) regularisation
for (s in 1:S) {
  prior_ += dirichlet_lpdf(B[s] | rep_vector(2, M));
}
