// priors for (very weak) regularisation
prior_ += normal_lpdf(beta_o_raw | 0, 5);

