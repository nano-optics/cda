#include "utils.h"
#include "cda.h"
#include "linear.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;


arma::colvec linear_extinction(const arma::mat& R, const arma::cx_mat& A, const double kn)
{
  const int N = R.n_rows;
  //constants
  const arma::cx_double i = arma::cx_double(0,1);
  
  // incident field
  const arma::cx_colvec SP="(1,0) (0,0) (0,0);", PP="(0,0) (1,0) (0,0);";
  arma::cx_colvec Eincident(3*N), P(3*N);
  const arma::colvec  khat="0 0 1;", kvec = kn*khat;
  const arma::mat kr = R * kvec;
  const arma::cx_mat expikr = exp(i*kr);
  arma::colvec res(2);

  // s polarisation
  Eincident = reshape(expikr * strans(SP), 3*N, 1, 1);
  // P = A * Eincident;
  P = solve(A, Eincident);
  res(0) = extinction(kn, P, Eincident) / N; 

  // p polarisation
  Eincident = reshape(expikr * strans(PP), 3*N, 1, 1);
  // P = A * Eincident;
  P = solve(A, Eincident);
  res(1) = extinction(kn, P, Eincident) / N;  
  
  return res ;
  } 
 

arma::mat linear_extinction_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				     const arma::mat& Euler, const int progress)
{
  
  int N = kn.n_elem, Nr = R.n_rows, ll;
  // Rcpp::Rcout << N << "\n";
  arma::mat res(N,2);
  arma::cx_mat beta(3,Nr);
  arma::colvec tmp(2);
  arma::cx_mat A(3*Nr,3*Nr);
  
  for(ll=0; ll<N; ll++){ // loop over kn   
    if(progress == 1)
      progress_bar(ll+1,N);
      beta = reshape(Beta.row(ll), 3, Nr, 1); 
      A = interaction_matrix(R, kn[ll], beta, Euler, 1);
      tmp = linear_extinction(R, A, kn[ll]);
      res(ll,0) = tmp(0);
      res(ll,1) = tmp(1);
  }
  if(progress == 1)
    Rcpp::Rcout <<"\n";
  
  return res ;
} 


// // calculation for multiple angles of incidence
// in progress, not functional
// R: positions
// A: interaction matrix
// invalpha: inverse polarisabilities
// kn: wavevector
// Angles: incident angles
 arma::mat dispersion(const arma::mat& R, const arma::cx_mat& A, const arma::cx_mat& invalpha, \
		      const double kn, const arma::mat& Angles, const arma::mat& Euler, const int invert)
   {
     const int N = R.n_rows, NAngles = Angles.n_rows;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::math::pi();
    arma::mat Rot(3,3);

    arma::cx_mat polar = diagonal_polarisability(invalpha, Euler);

    // incident field
    const arma::cx_colvec LPP="(1,0) (0,0) (0,0);", LPS="(0,0) (1,0) (0,0);";
    arma::cx_colvec ELPP(3), ELPS(3), Eincident(3*N), P(3*N);
    const arma::colvec  khat="1 0 0;", kvec = kn*khat; 
    double phi, psi,theta;
    arma::mat kr;
    arma::cx_mat expikr, B;

    if(invert == 1){
     B = pinv(A); /* inverting the interaction matrix 
     				 to solve AP=Eincident multiple times */
    }
    arma::mat res(NAngles, 4) ;  

    // begin calculation

    int ll=0; 
    for(ll=0; ll<NAngles; ll++){ // loop over angles

      phi = Angles(ll, 0), psi = Angles(ll, 1), theta = Angles(ll, 2);
     
      Rot = euler(phi, pi/2, psi); // theta is fixed
      ELPP =  trans(Rot) * LPP ;
      ELPS =  trans(Rot) * LPS ;
      kr = R * trans(Rot) * kvec;
      expikr = exp(i*kr);
      
      // P polarisation
      Eincident = reshape(expikr * strans(ELPP), 3*N, 1, 1);
      
      if(invert == 1){
	P = B * Eincident;
      } else {
	P = solve(A, Eincident);
      }

      res(ll,0) =  extinction(kn, P, Eincident); 
      res(ll,1) =  absorption(kn, P, polar); 
      
      // S polarisation
      Eincident = reshape(expikr * strans(ELPS), 3*N, 1, 1);

      if(invert == 1){
	P = B * Eincident;
      } else {
	P = solve(A, Eincident);
      }

      res(ll,2) =  extinction(kn, P, Eincident); 
      res(ll,3) = absorption(kn, P, polar); 
      
    } 
             
    return res ;
   } 

arma::cube dispersion_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				      const arma::mat& Euler, const arma::mat& Angles, \
				      const int invert, const int progress)
  {

    const int NAngles = Angles.n_rows;
    int N = kn.n_elem, Nr = R.n_rows, ll;
    arma::cube res(NAngles, 4, N);
    arma::cx_mat beta(3,Nr);
    arma::mat tmp(NAngles, 4);
    arma::cx_mat A(3*Nr,3*Nr), polar(3*Nr,3*Nr);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress == 1)
	progress_bar(ll+1,N);
      beta = reshape(Beta.row(ll), 3, Nr, 1); 
      A = interaction_matrix(R, kn[ll], beta, Euler, 1); // always full
      // polar = diagonal_polarisability(beta, Euler);
      // tmp = dispersion(R, A, polar, kn[ll], Angles);
      tmp = dispersion(R, A, beta, kn[ll], Angles, Euler, invert);

      // cout << A << "\n";
      res.slice(ll) = tmp; 
    }
    if(progress == 1)
      Rcpp::Rcout << "\n";

    return res ;
  } 



RCPP_MODULE(linear){
       using namespace Rcpp ;

       function( "linear_extinction_spectrum", &linear_extinction_spectrum, \
		 "Returns the extinction spectra for x and y polarisation at fixed incidence" ) ;
       function( "dispersion", &dispersion, \
       		 "Returns the abs and ext xsec for x and y polarisation at multiple angles of incidence" ) ;
       function( "dispersion_spectrum", &dispersion_spectrum,		\
       		 "Returns the abs and ext xsec for x and y polarisation at multiple angles of incidence" ) ;
}
