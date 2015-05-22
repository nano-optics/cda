#ifndef _cda_CG_H
#define _cda_CG_H

#include <RcppArmadillo.h>

arma::cx_colvec cg_solve(const arma::cx_mat& A, 
			const arma::cx_vec& y,	
			const arma::cx_vec& x0,		
			const double nmax, 	
			const double tol);
#endif
