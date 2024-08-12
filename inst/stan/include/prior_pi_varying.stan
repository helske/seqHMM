// priors for (very weak) regularisation
prior_ += normal_lpdf(alpha_i_raw | 0, 5);
prior_ += normal_lpdf(to_vector(beta_i_raw) | 0, 5);
