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


// given a 3xN matrix Beta of principal (inverse) polarizabilities, and a
// 3xN matrix Euler of rotation angles
// returns a 3Nx3N block-diagonal matrix of (inverse) polarizabilities
// which are the diagonal blocks of the interaction matrix
arma::cx_mat block_diagonal(const arma::cx_mat& Beta, const arma::mat& Euler) {
  
  const int N = Euler.n_rows;
  const arma::colvec phi = Euler.col(0), theta = Euler.col(1), psi = Euler.col(2);
  
  arma::mat Rot(3,3);
  arma::cx_mat polar(3*N,3*N);
  
  int ii=0;
  for(ii=0; ii<N; ii++){
    
    Rot = euler(phi[ii], theta[ii], psi[ii]);
    polar.submat(ii*3,ii*3,ii*3+2,ii*3+2) = inv(Rot) * diagmat(Beta.col(ii)) * Rot; 
    
  } // polar is done
  
  return polar;
}

// constructs the interaction matrix from a Nx3 matrix of positions R,
// a wavenumber kn, a 3xN matrix Beta of principal inverse polarizabilities, a
// 3xN matrix Euler of rotation angles, a flag full to use the full
// dipole field or only the static term
arma::cx_mat interaction_matrix(const arma::mat& R, const double kn,	\
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
  arma::cx_mat polar = block_diagonal(Beta, Euler);
  
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
// polarization P, 3Nx3N block diagonal matrix diagBeta of inverse polarizabilities
double absorption(const double kn, const arma::cx_colvec& P, const arma::cx_mat& diagBeta)
{
  const double c = 4*arma::math::pi()*kn*(imag(cdot(diagBeta * P, P)) -	\
					  kn*kn*kn* 2/3 * real(cdot(P, P))); 
  return c;
}



// // calculation for multiple angles of incidence
// in progress, not functional
// R: positions
// A: interaction matrix
// invalpha: inverse polarisabilities
// kn: wavevector
// Angles: incident angles
 arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, const arma::cx_mat& invalpha, \
		      const double kn, const arma::vec& Angle, const arma::uvec& Axis,  \
		      const arma::mat& Euler, \
		      const int polarisation,			\
		      const int invert)
   {
     const int N = R.n_rows, NAngles = Angle.n_elem;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::math::pi();
    arma::mat Rot(3,3);

    arma::cx_mat polar = block_diagonal(invalpha, Euler);

    // incident field
    arma::cx_colvec LPP, LPS;
    
    if(polarisation == 0){ // linear
     LPP="(1,0) (0,0) (0,0);", LPS="(0,0) (1,0) (0,0);";
      } else { // circular
     LPP="(0,1) (1,0) (0,0);", LPS="(1,0) (0,1) (0,0);";
    }

    const arma::colvec  khat="0 0 1;"; 
    const arma::colvec kvec = kn*khat; 
   

    arma::cx_colvec ELPP(3), ELPS(3), Eincident(3*N), P(3*N);
    double phi, psi,theta;
    arma::mat kr;
    arma::cx_mat expikr, B;

    if(invert == 1){
     B = pinv(A); /* inverting the interaction matrix 
     				 to solve AP=Eincident multiple times */
    }
    arma::mat res(NAngles, 6) ;  

    // begin calculation

    int ll=0; 
    for(ll=0; ll<NAngles; ll++){ // loop over angles

      Rot = axis_rotation(Angle(ll), Axis(ll)); // note: rotation along fixed axis only
      ELPP =  trans(Rot) * LPP ;
      ELPS =  trans(Rot) * LPS ;
      kr = R * trans(Rot) * kvec;
      expikr = exp(i*kr);
      
      // PL or RC polarisation
      Eincident = reshape(expikr * strans(ELPP), 3*N, 1, 1);
      
      if(invert == 1){
	P = B * Eincident;
      } else {
	P = solve(A, Eincident);
      }

      res(ll,0) =  extinction(kn, P, Eincident); 
      res(ll,1) =  absorption(kn, P, polar); 
      res(ll,2) =  res(ll,0) - res(ll,1); 
      
      // SL or LC polarisation
      Eincident = reshape(expikr * strans(ELPS), 3*N, 1, 1);

      if(invert == 1){
	P = B * Eincident;
      } else {
	P = solve(A, Eincident);
      }

      res(ll,3) =  extinction(kn, P, Eincident); 
      res(ll,4) = absorption(kn, P, polar); 
      res(ll,5) =  res(ll,3) - res(ll,4); 
      
    } 
             
    return res ;
   } 

arma::cube dispersion_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
			       const arma::mat& Euler, const arma::vec& Angle, const arma::uvec& Axis, \
			       const int polarisation,			\
			       const int invert, const int progress)
  {

    const int NAngles = Angle.n_elem;
    int N = kn.n_elem, Nr = R.n_rows, ll;
    arma::cube res(NAngles, 6, N);
    arma::cx_mat beta(3,Nr);
    arma::mat tmp(NAngles, 6);
    arma::cx_mat A(3*Nr,3*Nr), polar(3*Nr,3*Nr);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress == 1)
	progress_bar(ll+1,N);
      beta = reshape(Beta.row(ll), 3, Nr, 1); 
      A = interaction_matrix(R, kn[ll], beta, Euler, 1); // always full
      tmp = dispersion(R, A, beta, kn[ll], Angle, Axis, Euler, polarisation, invert);

      res.slice(ll) = tmp; 
    }
    if(progress == 1)
      Rcpp::Rcout << "\n";

    return res ;
  } 


RCPP_MODULE(cda){
       using namespace Rcpp ;

       function( "euler", &euler, "Constructs a 3x3 Euler rotation matrix" ) ;
       function( "extinction", &extinction, "Calculates the extinction cross-section" ) ;
       function( "absorption", &absorption, "Calculates the absorption cross-section" ) ;
       function( "interaction_matrix", &interaction_matrix, "Constructs the coupled-dipole interaction matrix" ) ;
       function( "dispersion_spectrum", &dispersion_spectrum,		\
       		 "Returns the absorption, scattering and extinction spectra for 2 polarisations at specified angles of incidence" ) ;
}
