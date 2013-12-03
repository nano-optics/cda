// 
// Multiple angles of incidence
// 
#include "utils.h"
#include "cda.h"
#include "dispersion.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;
using namespace arma;


//
// Angle-resolved cross-sections for multiple directions of incidence
//
// R is the Nx3 matrix of positions
// A is the 3Nx3N interaction matrix
// Adiag is the 3Nx3N block-diagonal part of the interaction matrix
// kn is the incident wavenumber (scalar)
// Angles is the Nangles vector of incident beam angles
// Axes is the Nangles vector of incident beam axes
// polarisation is an integer flag to switch between linear and circular polarisation
 arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, 
		      const arma::cx_mat& Adiag,			
		      const double kn, const arma::vec& Angles, 
		      const arma::ivec& Axes, 
		      const int polarisation)
   {
     const int N = R.n_rows, NAngles = Angles.n_elem;
     arma::mat Rot(3,3);


    // incident field
    const arma::colvec  khat="0 0 1;"; 
    const arma::colvec kvec = kn*khat; 
    arma::cx_colvec LPP, LPS;
    if(polarisation == 0){ // linear
      LPP="(1,0) (0,0) (0,0);", LPS="(0,0) (1,0) (0,0);";
    } else { // circular
      LPP="(1,0) (0,1) (0,0);", LPS="(0,1) (1,0) (0,0);";
      LPP = arma::datum::sqrt2/2 * LPP ;
      LPS = arma::datum::sqrt2/2 * LPS ;
    }

    arma::mat res(NAngles, 6) ;  
    arma::cx_mat Eincident(3*N,NAngles);
    arma::cx_mat P(3*N,NAngles);

    // first polarisation
    Eincident = multiple_incident_field(LPP, kvec, R, Axes, Angles);
    P = solve(A, Eincident);
    res.col(0) =  extinction(kn, P, Eincident); 
    res.col(1) = absorption(kn, P, Adiag); 
    res.col(2) = res.col(0) - res.col(1); 

    // second polarisation
    Eincident = multiple_incident_field(LPS, kvec, R, Axes, Angles);
    P = solve(A, Eincident);
    res.col(3) =  extinction(kn, P, Eincident); 
    res.col(4) = absorption(kn, P, Adiag); 
    res.col(5) = res.col(3) - res.col(4); 
             
    return res ;
   } 


//
// Angle-resolved spectra for linear or circular polarisations
//
// kn is the vector of incident wavenumbers
// Beta is the 3N vector of inverse polarisabilities
// R is the Nx3 matrix of positions
// Euler is the Nx3 matrix of particle rotation angles
// Angles is the Nangles vector of incident beam angles
// Axes is the Nangles vector of incident beam axes
// polarisation is an integer flag to switch between linear and circular polarisation
// progress is a logical flag to display progress bars
arma::cube dispersion_spectrum(const arma::colvec kn, 
			       const arma::cx_mat& Beta, const arma::mat& R, 
			       const arma::mat& Euler, const arma::vec& Angles,
			       const arma::ivec& Axes,			
			       const int polarisation, const bool progress)
  {

    const int NAngles = Angles.n_elem;
    int N = kn.n_elem, Nr = R.n_rows, ll;
    arma::cube results(NAngles, 6, N);
    arma::mat tmp(NAngles, 6);
    arma::cx_mat A(3*Nr,3*Nr), Adiag(3*Nr,3*Nr);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress)
	progress_bar(ll+1,N);
      A = interaction_matrix(R, kn(ll), Beta.col(ll), Euler, 1); // retarded
      Adiag = block_diagonal(Beta.col(ll), Euler);
      tmp = dispersion(R, A, Adiag, kn(ll), Angles, Axes, polarisation);
      results.slice(ll) = tmp; 
    }
    if(progress)
      Rcpp::Rcout << "\n";

    return results ;
  } 


RCPP_MODULE(dispersion){

       Rcpp::function( "dispersion", &dispersion,
	    "Angle-resolved cross-sections for multiple directions of incidence" ) ;
       Rcpp::function( "dispersion", &dispersion,
	    "Angle-resolved cross-sections for multiple directions of incidence" ) ;
       Rcpp::function( "dispersion_spectrum", &dispersion_spectrum,
	   "Angle-resolved spectra for linear or circular polarisations" ) ;
}
