//log-sum-exp trick
#include "logsumexp.h"

// No idea how I decided to implement logSumExp like this, but apparently this
// is now recommended in
// Pierre Blanchard, Desmond J Higham, Nicholas J Higham (2021). 
// Accurately computing the log-sum-exp and softmax functions, 
// IMA Journal of Numerical Analysis, 41, 4, 2311–2330
// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
  arma::uword maxi = x.index_max();
  double maxv = x(maxi);
  if (!std::isfinite(maxv)) {
    return maxv;
  }
  double cumsum = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}
