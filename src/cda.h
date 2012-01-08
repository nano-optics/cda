#ifndef _cda_CDA_H
#define _cda_CDA_H

#include <RcppArmadillo.h>

arma::cx_mat diagonal_polarisability(const arma::cx_mat& Alpha, const arma::mat& Euler);
arma::cx_mat interaction_matrix(const arma::mat& R, const double kn, const arma::cx_mat& invAlpha, const arma::mat& Euler, const int full);
double extinction(const double kn, const arma::cx_colvec& P, const arma::cx_colvec& Eincident);
double absorption(const double kn, const arma::cx_colvec& P, const arma::cx_mat& invpolar);

#endif
