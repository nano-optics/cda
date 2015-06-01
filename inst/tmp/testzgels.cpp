#define ARMA_EXTRA_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cx_mat main_call(const arma::cx_mat& A, const arma::cx_mat& b)
{
  const int n=A.n_rows;
  arma::cx_mat x;
  //x = solve(A, b);
  bool success, slow=true;
  success =  arma::auxlib::solve(x, A, b, slow);
  return x;
} 


/*** R
set.seed(123)
n <- 10
Asquare = matrix(rnorm(n^2)+0i, n)
Arect = matrix(rnorm(n*3)+0i, nrow = n)
b <- matrix(rnorm(n) + 0i, ncol=1)

x1 <- main_call(Asquare, b)
# ... gesv is called (LU decomposition)
x2 <- main_call(Arect, b)
# ... gels is called (least-square solution of over-determined system)
*/
