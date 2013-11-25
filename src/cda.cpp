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

//
// Calculate the incident field at each dipole location 
//
// for multiple (Euler) angles of incidence
// E0 is the normalised incident electric field before rotation
// k is the incident wavevector before rotation
// R is the Nx3 matrix of positions
// Angles is the Nx3 matrix of incident beam angles
arma::cx_mat incident_field(const arma::cx_colvec& E0, 
			    const arma::colvec& k, 
			    const arma::mat& R,
			    const arma::mat& Angles)
{
  const int Nangles = Angles.n_rows;
  const int N = R.n_rows;
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
  for(jj=0; jj<Nangles; jj++)
    {
      Rot = euler(Angles(jj,0), Angles(jj,1), Angles(jj,2));
      k_r = Rot.st() * k;
      E0_r = Rot.st() * E0;
      kR = R * k_r ;
      expikr = exp(i * kR);
      expikrrep = strans(vectorise(repmat(expikr, 1, 3), 1));
      E0rep = repmat(E0_r, N, 1);
      Ei.col(jj) = E0rep % expikrrep;
    }
  return(Ei);
}

//
// Incident field along multiple axes
//
// E0 is the normalised electric field
// k is the wavevector
// R is the Nx3 matrix of positions
// Axes is a vector of integer codes corresponding to x, y, z
// Angles is a vector of rotation angles around Axes
arma::cx_mat multiple_incident_field(const arma::cx_colvec& E0, 
			     const arma::colvec& k, 
			     const arma::mat& R,
			     const arma::ivec& Axes,
			     const arma::colvec& Angles)
{
  const int Nangles = Angles.n_elem;
  const int N = R.n_rows;
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
  for(jj=0; jj<Nangles; jj++)
    {
      Rot = axis_rotation(Angles(jj), Axes(jj));
      k_r = Rot.st() * k;
      E0_r = Rot.st() * E0;
      kR = R * k_r ;
      expikr = exp(i * kR);
      expikrrep = strans(vectorise(repmat(expikr, 1, 3), 1));
      E0rep = repmat(E0_r, N, 1);
      Ei.col(jj) = E0rep % expikrrep;
    }
  return(Ei);
}

//
// Construct the diagonal blocks of the interaction matrix
// (inverse polarisabilities) for use in C_abs
//
// Beta is the 3N vector of inverse polarisabilities
// Nx3 matrix Euler of rotation angles
// returns a 3Nx3N block-diagonal matrix of (inverse) polarizabilities
// which are the diagonal blocks of the interaction matrix
arma::cx_mat block_diagonal(const arma::cx_colvec& Beta, const arma::mat& Euler) {
  
  const int N = Euler.n_rows;
  arma::mat Rot(3,3);
  arma::cx_mat polar = arma::zeros<arma::cx_mat>( 3*N, 3*N );
  int ii=0;
  for(ii=0; ii<N; ii++){
    Rot = euler(Euler(ii,0), Euler(ii,1), Euler(ii,2));
    polar.submat(ii*3,ii*3,ii*3+2,ii*3+2) =  Rot.st() * \
      diagmat(Beta.subvec(ii*3, ii*3+2)) * Rot; 
    
  } // polar is done
  
  return polar;
}

//
// Construct the full interaction matrix
//
// R is the Nx3 matrix of positions
// kn is the incident wavenumber (scalar)
// Beta is the 3N vector of inverse polarisabilities
// Euler is the Nx3 matrix of rotation angles
// full is an logical flag to use full/static interaction term
arma::cx_mat interaction_matrix(const arma::mat& R, const double kn,
				const arma::cx_colvec& Beta, const arma::mat& Euler, 
				const bool full) {
  
  const int N = R.n_rows;
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
	      rk_to_rj = R.row(jj) - R.row(kk) ;
	      rjk = norm(rk_to_rj, 2);
	      rjkhat = rk_to_rj / rjk;
	      rjkrjk =  rjkhat.st() * rjkhat;
	      if(full) {
		Ajk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) + \
					      (i*kn*rjk - arma::cx_double(1,0)) / \
					      (rjk*rjk) * (3*rjkrjk - I3)) ;
	      } else {	      
		Ajk = (I3 - 3*rjkrjk)/ (rjk*rjk*rjk)  ;
	      }
	    } else { // diagonal blocks
	    Rot = euler(Euler(jj,0), Euler(jj,1), Euler(jj,2));
	    Ajk = Rot.st() * diagmat(Beta.subvec(jj*3, jj*3+2)) * Rot;
	      
	      }
	  // assign block 
	  A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk;
	}
    } // end loops
  
  return(A);
}


//
// Calculate the extinction cross-section for multiple incident angles
//
// wavenumber kn
// Nx3 polarisation P
// Nx3 incident field Eincident
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


//
// Calculate the absorption cross-section for multiple incident angles
//
// wavenumber kn
// Nx3 polarisation P
// 3Nx3N block diagonal matrix diagBeta of inverse polarizabilities
 arma::colvec absorption(const double kn, const arma::cx_mat& P, 
		  const arma::cx_mat& Adiag)
{
  int Nangles = P.n_cols, ii=0;
  arma::cx_mat Eexc = Adiag * P;
  arma::colvec results(Nangles);

  for (ii=0; ii<Nangles; ii++)
    {
      results(ii) = imag(cdot(Eexc.col(ii), P.col(ii))) -	\
			      kn*kn*kn* 2/3 * real(cdot(P.col(ii), P.col(ii)));
    }
  return  4*arma::datum::pi*kn*results;

}


RCPP_MODULE(cda){
       Rcpp::function( "euler", &euler, 
		       "Euler rotation matrix" ) ;
       Rcpp::function( "axis_rotation", &axis_rotation, 
		       "Rotation matrix about a cartesian axis" ) ;
       Rcpp::function( "extinction", &extinction, 
	 "Calculate the extinction cross-section for multiple incident angles" ) ;
       Rcpp::function( "absorption", &absorption, 
         "Calculate the absorption cross-section for multiple incident angles" ) ;
       Rcpp::function( "interaction_matrix", &interaction_matrix, 
		       "Construct the full interaction matrix" ) ;
       Rcpp::function( "block_diagonal", &block_diagonal, 
		       "Construct the diagonal blocks of the interaction matrix" ) ;
       Rcpp::function( "incident_field", &incident_field, 
		       "Calculate the incident field at each dipole location" ) ;
       Rcpp::function( "multiple_incident_field", &multiple_incident_field,
	    "Incident field along multiple axes" ) ;
}
