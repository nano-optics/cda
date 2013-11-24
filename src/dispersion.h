#ifndef _cda_DISPERSION_H
#define _cda_DISPERSION_H

#include <RcppArmadillo.h>

 arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, 
		      const arma::cx_mat& Adiag,			
		      const double kn, const arma::vec& Angles, 
		      const arma::ivec& Axes, 
		      const int polarisation);

#endif
