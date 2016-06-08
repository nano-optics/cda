#ifndef _AVERAGING_H
#define _AVERAGING_H

arma::rowvec cpp_oa(const arma::mat& R,
                    const arma::cx_mat& A,
                    const arma::cx_cube& AlphaBlocks,
                    const double kn,
                    const double medium,
                    const arma::mat& IncidenceNodes,
                    const arma::colvec& IncidenceWeights,
                    const arma::mat& ScatteringNodes,
                    const arma::colvec& ScatteringWeights,
                    const int inversion,
                    const int maxiter, const double tol);

#endif
