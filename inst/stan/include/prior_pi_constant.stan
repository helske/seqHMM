// priors for (very weak) regularisation
prior_ += dirichlet_lpdf(Pi | rep_vector(2, S));
