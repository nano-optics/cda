#ifndef _OOS_H
#define _OOS_H

arma::colvec cpp_iterate_field(const arma::mat& R,
                               const double kn, const double medium,
                               const arma::cx_mat& Ein,
                               const arma::cx_cube& AlphaBlocks,
                               arma::cx_mat& Eloc,
                               arma::cx_mat& P);

#endif
