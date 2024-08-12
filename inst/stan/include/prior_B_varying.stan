// priors for (very weak) regularisation
for(s in 1:S) {
  prior_ += normal_lpdf(to_vector(beta_o_raw[s]) | 0, 5);
}
