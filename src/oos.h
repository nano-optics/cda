#ifndef _OOS_H
#define _OOS_H

arma::colvec cpp_iterate_field(const arma::mat& R,
  const double kn, const double medium,
  const arma::cx_mat& Ein,
  const arma::cx_cube& AlphaBlocks,
  arma::cx_mat& Eloc,
  arma::cx_mat& P);

  arma::rowvec cpp_oa_oos(const arma::mat& R,
    const arma::cx_cube& AlphaBlocks,
    double kn, double medium,
    const arma::mat& IncidenceNodes,
    const arma::colvec& IncidenceWeights,
    const arma::mat& ScatteringNodes,
    const arma::colvec& ScatteringWeights,
    int maxiter, double tol);

#endif
