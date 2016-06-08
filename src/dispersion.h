#ifndef _DISPERSION_H
#define _DISPERSION_H

arma::mat cpp_dispersion(const arma::mat& R,
                         const arma::cx_mat& A,
                         const arma::cx_cube& AlphaBlocks,
                         const double kn, const double medium,
                         const arma::vec& Incidence,
                         const arma::ivec& Axes,
                         const arma::mat& ScatteringNodes,
                         const arma::colvec& ScatteringWeights,
                         const int polarisation,
                         const int inversion,
                         const int maxiter, const double tol);

#endif
