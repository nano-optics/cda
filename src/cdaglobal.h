#ifndef _CDAGLOBAL_H
#define _CDAGLOBAL_H

void cpp_interaction_matrix_update(const arma::mat& R, const double kn,
                                   const arma::cx_cube& AlphaBlocks,
                                   arma::cx_mat& A) ;

void cpp_propagator_update(const arma::mat& R, const double kn,
                                   const arma::cx_cube& AlphaBlocks,
                                   arma::cx_mat& G) ;

void cpp_polarization_update(const arma::cx_mat& E,
                             const arma::cx_cube& AlphaBlocks,
                             arma::cx_mat& P) ;

void cpp_alpha_blocks_update(const arma::cx_colvec& Alpha,
                             const arma::mat& Angles,
                             arma::cx_cube& AlphaBlocks) ;

#endif
