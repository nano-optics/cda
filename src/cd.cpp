#include "utils.h"
#include "cda.h"
#include "cd.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

// optional arguments in c++
// switch between c++ functions
// Toeplitz
// Gauss Legendre quadrature

// angular averaging with unit weights
arma::colvec averaging(const arma::mat& R, const arma::cx_mat& A, 
			const arma::cx_mat& Beta,		
			const double kn, 
			const arma::mat& Angles,
			const arma::colvec& Weights)
  {

    const int N = R.n_cols, Nangles = Angles.n_cols;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::datum::pi;
    arma::colvec res(4) ;   

    // incident field
    arma::cx_mat Eincident(3*N, Nangles), P(3*N, Nangles);
    const arma::colvec  khat="1 0 0;", kvec = kn*khat;
    arma::cx_colvec RCP="(0,0) (0,1) (1,0);", LCP="(0,0) (1,0) (0,1);";
    RCP = arma::datum::sqrt2/2 * RCP ;
    LCP = arma::datum::sqrt2/2 * LCP ;
    arma::colvec xsec(Nangles); // temporary storage of cross-sections

    // left polarisation
    Eincident = incident_field(LCP, kvec, R, Angles);
    P = solve(A, Eincident);
    xsec = extinction(kn, P, Eincident);
    res(0) = dot(xsec, Weights); 
    xsec = absorption(kn, P, Beta);
    res(1) = dot(xsec, Weights); 
      
    // right polarisation
    Eincident = incident_field(RCP, kvec, R, Angles);
    P = solve(A, Eincident);
    xsec = extinction(kn, P, Eincident);
    res(2) = dot(xsec, Weights); 
    xsec = absorption(kn, P, Beta);
    res(3) = dot(xsec, Weights); 

    return res ;
  } 

RCPP_MODULE(cd){
       Rcpp::function( "averaging", &averaging, \
		 "Calculates the orientation-averaged spectrum for absorption and extinction using numerical quadrature over two Euler angles" ) ;


}






