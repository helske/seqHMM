#include <algorithm>
#include <cmath>

double logSumExp(const double& x, const double& y) {
  double maxxy = std::max(x,y);
  double minxy = std::min(x,y);
  return maxxy + log1p(exp(minxy-maxxy));
}