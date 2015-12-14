//Armadillo version of R's rep(x, times)
#include <algorithm>
#include "seqHMM.h"
arma::vec reparma(const arma::vec& x, const arma::ivec& y) {
  arma::vec res(sum(y));
  int ind = 0;
  for (int i = 0; i < y.n_elem; ++i) {
    std::fill(res.begin() + ind, res.begin() + ind + y(i), x(i));
    ind += y(i);
  }
  return res;
}
