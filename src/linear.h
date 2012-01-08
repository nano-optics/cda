#ifndef _cda_LINEAR_H
#define _cda_LINEAR_H

#include <RcppArmadillo.h>

arma::colvec linear_extinction(const arma::mat& R, const arma::cx_mat& A, const double kn);
arma::mat linear_extinction_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				     const arma::mat& Euler, const int progress);

#endif
