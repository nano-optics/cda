#ifndef _CG_H
#define _CG_H

arma::cx_mat cpp_cg_solve(const arma::cx_mat& A,
                          const arma::cx_mat& y,
                          arma::cx_mat& x0,
                          const double nmax,
                          const double tol);

#endif
