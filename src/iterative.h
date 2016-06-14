#ifndef _OOS_H
#define _OOS_H

bool cpp_iterate_field(const arma::cx_mat& E0,
                      const arma::cx_mat& G,
                      const arma::cx_cube& AlphaBlocks,
                      const double kn,
                      const double tol,
                      const int maxiter,
                      arma::cx_mat& Eloc,
                      arma::cx_mat& P);

#endif
