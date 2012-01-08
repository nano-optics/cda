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
  // cout << N << "\n";
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
    cout<<"\n";
  
  return res ;
} 



RCPP_MODULE(linear){
       using namespace Rcpp ;

       function( "linear_extinction_spectrum", &linear_extinction_spectrum, \
		 "Returns the extinction spectra for x and y polarisation at fixed incidence" ) ;
    
}
