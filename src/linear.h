#ifndef _cda_LINEAR_H
#define _cda_LINEAR_H

#include <RcppArmadillo.h>

arma::colvec linear_extinction(const arma::mat& R, const arma::cx_mat& A, const double kn);
arma::mat linear_extinction_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				     const arma::mat& Euler, const int progress);

arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, const arma::cx_mat& invalpha, \
		     const double kn, const arma::mat& Angles,  const arma::mat& Euler, const int invert);
#endif
