#ifndef _cda_CG_H
#define _cda_CG_H

#include <RcppArmadillo.h>

arma::cx_mat cg_solve(const arma::cx_mat& A, 
			 const arma::cx_mat& y,	
			 arma::cx_mat& x0,		
			 const double nmax, 	
			 const double tol);
#endif
