        for(s in 1:S) {
          log_B[s, 1:M] = log_softmax(beta_o[d, s] * X_o[t, i])';
        }
