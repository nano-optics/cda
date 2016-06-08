#ifndef _CDA_H
#define _CDA_H

arma::cx_mat cpp_interaction_matrix(const arma::mat& R, const double kn,
                                    const arma::cx_cube& AlphaBlocks);

arma::cx_mat cpp_propagator(const arma::mat& R, const double kn,
                            const arma::cx_cube& AlphaBlocks);

arma::cx_cube cpp_alpha_blocks(const arma::cx_colvec& Alpha,
                               const arma::mat& Angles);

arma::cx_mat cpp_polarization(const arma::cx_mat& E,
                              const arma::cx_cube& AlphaBlocks);


#endif
