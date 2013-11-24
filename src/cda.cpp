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


// Beta is the 3N vector of inverse polarisabilities
// 3xN matrix Euler of rotation angles
// returns a 3Nx3N block-diagonal matrix of (inverse) polarizabilities
// which are the diagonal blocks of the interaction matrix
arma::cx_mat block_diagonal(const arma::cx_colvec& Beta, const arma::mat& Euler) {
  
  const int N = Euler.n_cols;
  arma::mat Rot(3,3);
  arma::cx_mat polar = arma::zeros<arma::cx_mat>( 3*N, 3*N );
  int ii=0;
  for(ii=0; ii<N; ii++){
    Rot = euler(Euler(0,ii), Euler(1,ii), Euler(2,ii));
    polar.submat(ii*3,ii*3,ii*3+2,ii*3+2) =  Rot.st() * \
      diagmat(Beta.subvec(ii*3, ii*3+2)) * Rot; 
    
  } // polar is done
  
  return polar;
}

// constructs the interaction matrix 
// R is the 3xN matrix of positions
// kn is the incident wavenumber (scalar)
// Beta is the 3N vector of inverse polarisabilities
// Euler is the 3xN matrix of rotation angles
// full is an integer flag to use full/static interaction term
arma::cx_mat interaction_matrix(const arma::mat& R, const double kn,
				const arma::cx_colvec& Beta, const arma::mat& Euler, 
				const int full) {
  
  const int N = R.n_cols;
  arma::cx_mat A = arma::zeros<arma::cx_mat>( 3*N, 3*N );

  // constants
  const arma::cx_double i = arma::cx_double(0,1);
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );

  // temporary vars
  arma::mat Rot(3,3);
  int jj=0, kk=0;
  arma::mat rk_to_rj = arma::mat(1,3), rjkhat = arma::mat(1,3) , 
            rjkrjk = arma::mat(3,3);
  
  double rjk;
  arma::cx_mat Ajk = arma::cx_mat(3,3);
  
  // nested for loop over dipole locations
  for(jj=0; jj<N; jj++)
    {
      for(kk=0; kk<N; kk++)
	{
	  if(jj!=kk)
	    {
	      rk_to_rj = R.col(jj) - R.col(kk) ;
	      rjk = norm(rk_to_rj,2);
	      rjkhat = rk_to_rj / rjk;
	      rjkrjk =  rjkhat * rjkhat.st();
	      if(full == 1) {
		Ajk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) + \
					      (i*kn*rjk - arma::cx_double(1,0)) / \
					      (rjk*rjk) * (3*rjkrjk - I3)) ;
	      } else {	      
		Ajk = (I3 - 3*rjkrjk)/ (rjk*rjk*rjk)  ;
	      }
	    } else { // diagonal blocks
	    Rot = euler(Euler(0,jj), Euler(1,jj), Euler(2,jj));
	    Ajk = Rot.st() * diagmat(Beta.subvec(jj*3, jj*3+2)) * Rot;
	      
	      }
	  // assign block 
	  A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk;
	}
    } // end loops
  
  return(A);
}


// calculate the extinction cross section given wavenumber kn, Nx3
// polarization P, Nx3 incident field Eincident
arma::colvec extinction(const double kn, const arma::cx_mat& P, 
			 const arma::cx_mat& Eincident)
{
  int Nangles = P.n_cols, ii=0;
  arma::colvec results(Nangles);

  for (ii=0; ii<Nangles; ii++)
    {
      results(ii) = imag(cdot(Eincident.col(ii), P.col(ii)));
    }
  return  4*arma::datum::pi*kn*results;
}


// calculate the absorption cross section given wavenumber kn, Nx3
// polarization P, 3Nx3N block diagonal matrix diagBeta of inverse polarizabilities
 arma::colvec absorption(const double kn, const arma::cx_mat& P, 
		  const arma::cx_mat& diagBeta)
{
  int Nangles = P.n_cols, ii=0;
  arma::cx_mat Eexc = diagBeta * P;
  arma::colvec results(Nangles);

  for (ii=0; ii<Nangles; ii++)
    {
      results(ii) = imag(cdot(Eexc.col(ii), P.col(ii))) -	\
			      kn*kn*kn* 2/3 * real(cdot(P.col(ii), P.col(ii)));
    }
  return  4*arma::datum::pi*kn*results;

}

// calculate the incident field at each dipole location 
// for multiple angles of incidence
arma::cx_mat incident_field(const arma::cx_colvec& E0, 
			    const arma::colvec& k, 
			    const arma::mat& R,
			    const arma::mat& Angles)
{
  const int Nangles = Angles.n_cols;
  const int N = R.n_cols;
  const arma::cx_double i = arma::cx_double(0,1);
  arma::mat Rot(3,3);
  arma::cx_mat Ei = arma::cx_mat(3*N,Nangles);
  arma::cx_colvec E0_r(3);
  arma::cx_colvec expikr(N);
  arma::colvec k_r(3);
  arma::colvec kR(N);
  arma::cx_colvec expikrrep(3*N);
  arma::cx_colvec E0rep(3*N);
  int jj=0;
  arma::mat  Rt = R.st();
  for(jj=0; jj<Nangles; jj++)
    {
      Rot = euler(Angles(0, jj), Angles(1, jj), Angles(2, jj));
      k_r = Rot.st() * k;
      E0_r = Rot.st() * E0;
      kR = Rt * k_r ;
      expikr = exp(i * kR);
      expikrrep = strans(vectorise(repmat(expikr, 1, 3), 1));
      E0rep = repmat(E0_r, N, 1);
      Ei.col(jj) = E0rep % expikrrep;
    }
  return(Ei);
}

RCPP_MODULE(cda){
       Rcpp::function( "euler", &euler, "Constructs a 3x3 Euler rotation matrix" ) ;
       Rcpp::function( "extinction", &extinction, "Calculates the extinction cross-section" ) ;
       Rcpp::function( "absorption", &absorption, "Calculates the absorption cross-section" ) ;
       Rcpp::function( "interaction_matrix", &interaction_matrix, "Constructs the coupled-dipole interaction matrix" ) ;
       Rcpp::function( "block_diagonal", &block_diagonal, "Diagonal part of the interaction matrix" ) ;
Rcpp::function( "incident_field", &incident_field, "Incident field at multiple angles" ) ;

}
