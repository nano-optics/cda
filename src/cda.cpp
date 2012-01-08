// 
// main functions for coupled dipole equations
// 
#include "utils.h"
#include "cda.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;


// given a 3xN matrix invAlpha of principal (inverse) polarizabilities, and a
// 3xN matrix Euler of rotation angles
// returns a 3Nx3N block-diagonal matrix of (inverse) polarizabilities
// which are the diagonal terms of the interaction matrix
arma::cx_mat diagonal_polarisability(const arma::cx_mat& Alpha, const arma::mat& Euler) {
  
  const int N = Euler.n_rows;
  const arma::colvec phi = Euler.col(0), theta = Euler.col(1), psi = Euler.col(2);
  
  arma::mat Rot(3,3);
  arma::cx_mat polar(3*N,3*N);
  
  int ii=0;
  for(ii=0; ii<N; ii++){
    
    Rot = euler(phi[ii], theta[ii], psi[ii]);
    polar.submat(ii*3,ii*3,ii*3+2,ii*3+2) = inv(Rot) * diagmat(Alpha.col(ii)) * Rot; 
    
  } // polar is done
  
  return polar;
}

// constructs the interaction matrix from a Nx3 matrix of positions R,
// a wavenumber kn, a 3xN matrix B of principal inverse polarizabilities, a
// 3xN matrix Euler of rotation angles, a flag full to use the full
// dipole field or only the static term
arma::cx_mat interaction_matrix(const arma::mat& R, const double kn,	\
				const arma::cx_mat& invAlpha, const arma::mat& Euler, 
				const int full) {
  
  const int N = R.n_rows;
  // temporary vars
  arma::mat Rot(3,3);
  const arma::cx_double i = arma::cx_double(0,1);
  arma::cx_mat A = arma::zeros<arma::cx_mat>( 3*N, 3*N );
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );
  
  int jj=0, kk=0;
  arma::mat rk_to_rj = arma::mat(1,3), rjkhat = arma::mat(1,3) , rjkrjk = arma::mat(3,3);
  double rjk;
  arma::cx_mat Ajk = arma::cx_mat(3,3);
      
  // nested for loop over dipole locations
  for(jj=0; jj<N; jj++)
    {
      for(kk=0; kk<N; kk++)
	{
	  if(jj!=kk){
	    
    	    rk_to_rj = R.row(jj) - R.row(kk) ;
    	    rjk = norm(rk_to_rj,2);
    	    rjkhat = rk_to_rj / rjk;
    	    rjkrjk = trans(rjkhat) * rjkhat;
	    if(full == 1) {
	      Ajk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) + (i*kn*rjk - arma::cx_double(1,0)) / (rjk*rjk) * (3*rjkrjk - I3)) ;
	    } else {	      
	      Ajk = (I3 - 3*rjkrjk)/ (rjk*rjk*rjk)  ;
	    }
	    // assign block 
	    A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk;
	  }
	}
    } // end loops
  
  // diagonal blocks
  arma::cx_mat polar = diagonal_polarisability(invAlpha, Euler);
  
  A = A + polar;
  // return inv(A); 
  return(A);
}

// calculate the extinction cross section given wavenumber kn, Nx3
// polarization P, Nx3 incident field Eincident
double extinction(const double kn, const arma::cx_colvec& P, const arma::cx_colvec& Eincident)
{
  const double c = 4*arma::math::pi()*kn * imag(cdot(Eincident, P)); 
  return c;
}

// calculate the absorption cross section given wavenumber kn, Nx3
// polarization P, 3Nx3N block diagonal matrix invalpha of inverse polarizabilities
double absorption(const double kn, const arma::cx_colvec& P, const arma::cx_mat& invpolar)
{
  const double c = 4*arma::math::pi()*kn*(imag(cdot(invpolar * P, P)) -	\
					  kn*kn*kn* 2/3 * real(cdot(P, P))); 
  return c;
}


RCPP_MODULE(cda){
       using namespace Rcpp ;

       function( "euler", &euler, "Constructs a 3x3 Euler rotation matrix" ) ;
       function( "extinction", &extinction, "Calculates the extinction cross-section" ) ;
       function( "absorption", &absorption, "Calculates the absorption cross-section" ) ;
       function( "interaction_matrix", &interaction_matrix, "Constructs the coupled-dipole interaction matrix" ) ;

}
