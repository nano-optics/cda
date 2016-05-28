// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include "cg.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;

//
// Solves a linear system iteratively via conjugate-gradient
//
// A: complex matrix
// y:  right-hand side
// x0:  initial guess
// nmax:  max number of iterations
// tol:  tolerance
// [[Rcpp::export]]
arma::cx_mat cpp_cg_solve(const arma::cx_mat& A,
			 const arma::cx_mat& y,
			 arma::cx_mat& x0,
			 const double nmax,
			 const double tol)
{

  int n = 0, ii = 0; // counter
  const int nr = y.n_cols; // number of right-hand side cols
  double rel_error = 1e4; // large value
  const int N = y.n_rows;

  // temporary variables
  double alpha_i, beta_i;
  const arma::cx_mat B = A.t(); // Hermitian transpose.
                                // Note: diag blocks not symmetric, so no easier way
  arma::cx_colvec y_i(N), z_i(N), g_i(N), p_i(N), w_i(N), v_i(N), x_i(N);
  arma::cx_colvec g_ip1(N), p_ip1(N), w_ip1(N), v_ip1(N), x_ip1(N);

  // loop over RHS cols
  for (ii=0; ii < nr; ii++){

    // initial step
    y_i = y.col(ii);
    z_i = B * y_i;
    x_i = x0.col(ii);
    g_i = z_i - B * A * x_i;
    p_i = g_i;
    w_i = A * x_i;
    v_i = A * p_i;

    n = 0; rel_error = 1e4;
    while ((n < nmax) && (rel_error > tol)){

      alpha_i = real(cdot(g_i, g_i) / cdot(v_i, v_i));
      x_ip1 = x_i + alpha_i * p_i;
      w_ip1 = w_i + alpha_i * v_i;
      g_ip1 = z_i - B * w_ip1;
      beta_i = real(cdot(g_ip1,g_ip1) / cdot(g_i,g_i));
      // updating
      p_i = g_ip1 + beta_i * p_i;
      v_i = A * p_i;
      rel_error = norm(x_ip1 - x_i) / norm(x_ip1) ;
      x_i = x_ip1;
      if(rel_error < tol) break; // skip cost of next line
      w_i = A * x_i;
      g_i = g_ip1;
      n++;
    }; // end while
    x0.col(ii) = x_i;
  }; //end loop

  return x0;
}
