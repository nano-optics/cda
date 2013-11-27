// 
// main functions for coupled dipole equations
// 
#include "utils.h"
#include "cda.h"
#include "internals.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

// calculate the extinction cross section given wavenumber kn, Nx3
// polarization P, Nx3 incident field Eincident
double extinction_single(const double kn, const arma::cx_mat& P, 
		  const arma::cx_mat& Eincident)
{
  const double c = 4*arma::math::pi()*kn * \
    imag(cdot(vectorise(Eincident), vectorise(P))) / P.n_cols; 
  return c;
}

// calculate the absorption cross section given wavenumber kn, Nx3
// polarization P, 3Nx3N block diagonal matrix diagBeta of inverse polarizabilities
double absorption_single(const double kn, const arma::cx_mat& P, 
		  const arma::cx_mat& diagBeta)
{
  arma::cx_colvec Pvec = vectorise(P, 0); 
  arma::cx_colvec Evec=vectorise(diagBeta * P, 0);
  const double c = 4*arma::math::pi()*kn*(as_scalar(imag(Evec.t() * Pvec)) - \
						    kn*kn*kn* 2/3 * \
						    real(cdot(Pvec, Pvec))); 
  return c/P.n_cols;
}

// constructs the interaction matrix from a Nx3 matrix of positions R,
// a complex wavevector kn, a Nx3 matrix Beta of principal inverse polarizabilities, a
// Nx3 matrix Euler of rotation angles, a flag full to use the full
// dipole field or only the static term
double eigenmodes(const arma::mat& R, const arma::cx_mat& knc,	\
				const arma::cx_mat& Beta, const arma::mat& Euler, 
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
  arma::cx_double kn = as_scalar(knc);
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
	      Ajk = exp(i*kn*rjk) / rjk * \
		(kn*kn*(rjkrjk - I3) +	  \
		 (i*kn*rjk - arma::cx_double(1,0)) / (rjk*rjk) * (3*rjkrjk - I3)) ;
	    } else {	      
	      Ajk = (I3 - 3*rjkrjk)/ (rjk*rjk*rjk)  ;
	    }
	    // assign block 
	    A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk;
	  }
	}
    } // end loops
  
  // diagonal blocks
  arma::cx_mat beta = block_diagonal(Beta, Euler);
  
  A = A + beta;
  // return inv(A); 
  return(abs(det(A)));
}

// constructs the green tensor matrix from a Ndipx3 matrix of dipole positions R,
// evaluated at Nprobex3 probe locations
// the wavevector kn can be complex
arma::cx_mat green_matrix(const arma::mat& Rdip, \
			  const arma::mat& Rprobe, 
			  const arma::cx_mat& kn) 
{
  
  const int Ndip = Rdip.n_rows, Nprobe = Rprobe.n_rows;
  // const arma::cx_vec kn = as<arma::cx_vec>(_kn);
  // temporary vars
  arma::mat Rot(3,3);
  const arma::cx_double i = arma::cx_double(0,1);
  // Nprobe rows Ndip cols
  arma::cx_mat A = arma::zeros<arma::cx_mat>( 3*Nprobe, 3*Ndip );
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );
  
  int jj=0, kk=0;
  arma::mat rk_to_rj = arma::mat(1,3), rjkhat = arma::mat(1,3) , rjkrjk = arma::mat(3,3);
  double rjk;
  arma::cx_mat Ajk = arma::cx_mat(3,3);
      
  // nested for loop over probe and dipole locations
  // assumes clean probe locations, distinct from dipoles

  for(jj=0; jj<Nprobe; jj++) // jj is for the probe
    {
      for(kk=0; kk<Ndip; kk++) // kk is for the dipoles
  	{

    	    rk_to_rj = Rprobe.row(jj) - Rdip.row(kk) ;
    	    rjk = norm(rk_to_rj, 2);
    	    rjkhat = rk_to_rj / rjk;
    	    rjkrjk = trans(rjkhat) * rjkhat;
	  
  	    Ajk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) + \
  					  (i*kn*rjk - arma::cx_double(1,0)) / \
  					  (rjk*rjk) * (3*rjkrjk - I3)) ;
	  
  	    // assign block 
  	    A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk;
  	}
	
    } // end loops
  
  return(A);
}


arma::cx_colvec polarization(const arma::mat& R, const arma::cx_mat& Beta, \
			     const arma::mat& Euler,	\
			     const double kn)
{
  const int Nr = R.n_rows;

  //constants
  const arma::cx_double i = arma::cx_double(0,1);
  
  // inverse polarizability
  arma::cx_mat beta(Nr,3);
  arma::cx_mat A(3*Nr,3*Nr);
  beta = reshape(Beta, Nr, 3, 0); 

  // incident field
  const arma::cx_colvec PP="(1,0) (0,0) (0,0);";
  const arma::colvec  khat="0 0 1;", kvec = kn*khat;
  const arma::mat kr = R * kvec;
  const arma::cx_mat expikr = exp(i*kr);

  A = interaction_matrix(R, kn, beta, Euler, 1);

  // p polarisation
  arma::cx_colvec Eincident(3*Nr), P(3*Nr);
  Eincident = reshape(expikr * strans(PP), 3*Nr, 1, 1);

  P = solve(A, Eincident);
 
  return P ;
  } 
 

// calculate the scattered field from a collection of dipoles P at Rdip,
// evaluated at Rprobe
arma::cx_mat  scattered_field(const arma::mat& Rdip, \
			      const arma::mat& Rprobe, const arma::cx_mat& kn, \
			      const arma::cx_colvec& P)
{
  // arma::cx_mat G = arma::zeros<arma::cx_mat>( 3*Nprobe, 3*Ndip );
  const arma::cx_mat G = green_matrix(Rdip, Rprobe, kn);
  const arma::cx_mat Esca = G * P;
  return Esca;
}

RCPP_MODULE(internals){

       Rcpp::function( "eigenmodes", &eigenmodes, "Computes the determinant of the interaction matrix" ) ;
       Rcpp::function( "green_matrix", &green_matrix, "Constructs the green matrix" ) ;
       Rcpp::function( "scattered_field", &scattered_field, "Calculates the scattered field" ) ;
       Rcpp::function( "polarization", &polarization, "Calculates the dipole moments" ) ;

}
