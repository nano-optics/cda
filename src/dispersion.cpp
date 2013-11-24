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


// calculate the incident field at each dipole location 
// for multiple axis-angles of incidence
arma::cx_mat multiple_incident_field(const arma::cx_colvec& E0, 
			     const arma::colvec& k, 
			     const arma::mat& R,
			     const arma::ivec& Axes,
			     const arma::colvec& Angles)
{
  const int Nangles = Angles.n_elem;
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
      Rot = axis_rotation(Angles(jj), Axes(jj));
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



// // calculation for multiple angles of incidence
// R: positions
// A: interaction matrix
// invalpha: inverse polarisabilities
// kn: wavevector
// Angles: incident angles
 arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, 
		      const arma::cx_mat& Adiag,			
		      const double kn, const arma::vec& Angles, 
		      const arma::ivec& Axes, 
		      const int polarisation)
   {
     const int N = R.n_cols, NAngles = Angles.n_elem;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::datum::pi;
    arma::mat Rot(3,3);


    // incident field
    const arma::colvec  khat="1 0 0;"; 
    const arma::colvec kvec = kn*khat; 
    arma::cx_colvec LPP, LPS;
    
    if(polarisation == 0){ // linear
     LPP="(0,0) (1,0) (0,0);", LPS="(0,0) (0,0) (1,0);";
      } else { // circular
      LPP="(0,0) (1,0) (0,1);", LPS="(0,0) (0,1) (1,0);";
      LPP = arma::datum::sqrt2/2 * LPP ;
      LPS = arma::datum::sqrt2/2 * LPS ;
    }

    arma::mat res(NAngles, 4) ;  
    arma::cx_mat Eincident(3*N,NAngles);
    arma::cx_mat P(3*N,NAngles);

    // first polarisation
    Eincident = multiple_incident_field(LPP, kvec, R, Axes, Angles);
    P = solve(A, Eincident);
    res.col(0) =  extinction(kn, P, Eincident); 
    res.col(1) = absorption(kn, P, Adiag); 

    // second polarisation
    Eincident = multiple_incident_field(LPS, kvec, R, Axes, Angles);
    P = solve(A, Eincident);
    res.col(2) =  extinction(kn, P, Eincident); 
    res.col(3) = absorption(kn, P, Adiag); 
             
    return res ;
   } 

RCPP_MODULE(dispersion){
       Rcpp::function( "dispersion", &dispersion,
		       "Returns the absorption, scattering and extinction spectra for 2 polarisations at specified angles of incidence" ) ;
}
