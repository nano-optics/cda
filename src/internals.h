#ifndef _cda_INTERNALS_H
#define _cda_INTERNALS_H

#include <RcppArmadillo.h>

arma::cx_mat green_matrix(const arma::mat& Rdip, \
			  const arma::mat& Rprobe, \
			  const arma::cx_mat& kn);

arma::cx_colvec polarization(const arma::mat& R, const arma::cx_mat& Beta, \
			  const arma::mat& Euler,			\
			  const double kn);
arma::cx_mat  scattered_field(const arma::mat& Rdip, \
			      const arma::mat& Rprobe, const arma::cx_mat& kn, \
			      const arma::cx_colvec& P);
#endif
