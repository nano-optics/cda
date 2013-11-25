#ifndef _cda_CDA_H
#define _cda_CDA_H

#include <RcppArmadillo.h>

arma::cx_mat block_diagonal(const arma::cx_colvec& Beta, const arma::mat& Euler);

arma::cx_mat interaction_matrix(const arma::mat& R, const double kn, 
				const arma::cx_colvec& Beta, const arma::mat& Euler,
				const bool full);

arma::cx_mat incident_field(const arma::cx_colvec& E0, 
			    const arma::colvec& k, 
			    const arma::mat& R,
			    const arma::mat& Angles);

arma::cx_mat multiple_incident_field(const arma::cx_colvec& E0, 
				     const arma::colvec& k, 
				     const arma::mat& R,
				     const arma::ivec& Axes,
				     const arma::colvec& Angles);

arma::colvec extinction(const double kn, const arma::cx_mat& P, 
				 const arma::cx_mat& Eincident);

arma::colvec absorption(const double kn, const arma::cx_mat& P, 
				  const arma::cx_mat& Adiag);


#endif
