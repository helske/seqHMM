// priors for (very weak) regularisation
for(s in 1:S) {
  prior_ += dirichlet_lpdf(A[s] | rep_vector(2, S));
}
