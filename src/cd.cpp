#include "utils.h"
#include "cda.h"
#include "cd.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;

using namespace std;

// angular averaging using Gauss-Legendre quadrature
arma::colvec averaging(const arma::mat& R, const arma::cx_mat& A, const arma::cx_mat& invalpha, \
			const double kn, const arma::mat& QuadPhi, const arma::mat& QuadPsi)
  {
    const int N = R.n_rows, NqPhi = QuadPhi.n_rows, NqPsi = QuadPsi.n_rows;
    const arma::colvec nodes1 = QuadPhi.col(0), weights1 = QuadPhi.col(1), \
      nodes2 = QuadPsi.col(0), weights2 = QuadPsi.col(1);;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::math::pi();
    arma::mat Rot(3,3);

    // incident field
    const arma::cx_colvec RCP="(0,0) (0,1) (1,0);", LCP="(0,0) (1,0) (0,1);";
    arma::cx_colvec ERCP(3), ELCP(3), Eincident(3*N), P(3*N);
    const arma::colvec  khat="1 0 0;", kvec = kn*khat;
    arma::mat kr;
    arma::cx_mat expikr;
    double left=0, right=0, left2=0, right2=0,phi=0,psi=0;
    
    // begin quadrature
    const double b2=pi/2, a2=2*pi, b1=-pi/2, a1=0; // a is for phi, b for psi
    const  double  C2 = (b2 - b1) / 2, D2 = (b2 + b1) / 2, \
      C1 = (a2 - a1) / 2,  D1 = (a2 + a1) / 2;
    
    const arma::colvec y = nodes2*C2 + D2; // for psi
    const arma::colvec x = nodes1*C1 + D1; // for phi
    
    double  tmpleft=0, tmpright=0, tmpleft2=0, tmpright2=0, factor, factor2;
    int ll=0,mm=0; // Nq quadrature points

    arma::cx_mat B = pinv(A); /* inverting the interaction matrix 
				 to solve AP=Eincident multiple times */

    for(ll=0; ll<NqPhi; ll++){ // loop over phi
      phi = x[ll];
      tmpleft=0; tmpright=0; tmpleft2=0; tmpright2=0;
	for(mm=0; mm<NqPsi; mm++){ // loop over psi
	  psi = y[mm];
	  Rot = euler(phi, pi/2, psi); // theta doesn't vary
	  ELCP = sqrt(2.0)/2 * trans(Rot) * LCP ;
	  ERCP = sqrt(2.0)/2 * trans(Rot) * RCP ;
	  kr = R * trans(Rot) * kvec;
	  expikr = exp(i*kr);
	  
	  factor = C1 * weights2[mm] * cos(psi) ;
	  // left polarisation
	  Eincident = reshape(expikr * strans(ELCP), 3*N, 1, 1);
	  P = B * Eincident;
	  // P = solve(A, Eincident);// too slow, invert A before loop
	  tmpleft +=  factor * extinction(kn, P, Eincident); 
	  tmpleft2 +=  factor * absorption(kn, P, invalpha); 
	  
	  // right polarisation
	  Eincident = reshape(expikr * strans(ERCP), 3*N, 1, 1);
	  P = B * Eincident; 
	  // P = solve(A, Eincident); // too slow, invert A before loop
	  tmpright +=  factor * extinction(kn, P, Eincident); 
	  tmpright2 +=  factor * absorption(kn, P, invalpha); 
	  
	}
	factor2 = C2 * weights1[ll];
	left  += factor2 * tmpleft;
	right += factor2 * tmpright;
	left2  += factor2 * tmpleft2;
	right2 += factor2 * tmpright2;
      } 

    arma::colvec res(4) ;               
                          
    res(0) = left / (4*pi); //ext L
    res(1) = right / (4*pi);//ext R
    res(2) = left2 / (4*pi); //abs L
    res(3) = right2 / (4*pi); //abs R
    return res ;
  } 

// angular averaging using Quasi Monte Carlo integration
arma::colvec averaging2(const arma::mat& R, const arma::cx_mat& A, const arma::cx_mat& invalpha, \
			const double kn, const arma::mat& QMC)
  {
    const int N = R.n_rows, NQMC = QMC.n_rows;
    //constants
    const arma::cx_double i = arma::cx_double(0,1);
    const double pi = arma::math::pi();
    arma::mat Rot(3,3);

    // incident field
    const arma::cx_colvec RCP="(0,0) (0,1) (1,0);", LCP="(0,0) (1,0) (0,1);";
    arma::cx_colvec ERCP(3), ELCP(3), Eincident(3*N), P(3*N);
    const arma::colvec  khat="1 0 0;", kvec = kn*khat;
    arma::mat kr;
    arma::cx_mat expikr;
    double left=0, right=0, left2=0, right2=0,phi=0,psi=0;
    
    arma::cx_mat B = pinv(A); /* inverting the interaction matrix 
				 to solve AP=Eincident multiple times */

    // begin integration
    
    int ll=0; 
    for(ll=0; ll<NQMC; ll++){ // loop over integration points
      phi = QMC(ll, 1)*2*pi, psi = asin(2*QMC(ll, 0) - 1);
     
      Rot = euler(phi, pi/2, psi); // theta doesn't vary
      ELCP = sqrt(2.0)/2 * trans(Rot) * LCP ;
      ERCP = sqrt(2.0)/2 * trans(Rot) * RCP ;
      kr = R * trans(Rot) * kvec;
      expikr = exp(i*kr);
      
      // left polarisation
      Eincident = reshape(expikr * strans(ELCP), 3*N, 1, 1);
      P = B * Eincident;
      // P = solve(A, Eincident);// too slow, invert A before loop
      left +=  extinction(kn, P, Eincident); 
      left2 +=  absorption(kn, P, invalpha); 
      
      // right polarisation
      Eincident = reshape(expikr * strans(ERCP), 3*N, 1, 1);
      P = B * Eincident; 
      // P = solve(A, Eincident); // too slow, invert A before loop
      right +=  extinction(kn, P, Eincident); 
      right2 += absorption(kn, P, invalpha); 
      
    } 

    arma::colvec res(4) ;               
    double normalisation;
    normalisation = NQMC;
    res(0) = left / normalisation; //ext L
    res(1) = right / normalisation;//ext R
    res(2) = left2 / normalisation; //abs L
    res(3) = right2 / normalisation; //abs R
    return res ;
  } 

arma::mat circular_dichroism_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				      const arma::mat& Euler, const arma::mat& QuadPhi, const arma::mat& QuadPsi, \
				      const int full, const int progress)
  {

    int N = kn.n_elem, Nr = R.n_rows, ll;
    // cout << N << "\n";
    arma::mat res(N,4);
    arma::cx_mat beta(3,Nr);
    arma::colvec tmp(4);
    arma::cx_mat A(3*Nr,3*Nr), polar(3*Nr,3*Nr);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress == 1)
	progress_bar(ll+1,N);
      beta = reshape(Beta.row(ll), 3, Nr, 1); 
      A = interaction_matrix(R, kn[ll], beta, Euler, full);
      polar = diagonal_polarisability(beta, Euler);
      tmp = averaging(R, A, polar, kn[ll], QuadPhi, QuadPsi);

      res(ll,0) = 0.5*(tmp(0) + tmp(1)); // extinction 
      res(ll,1) = 0.5*(tmp(2) + tmp(3)); // absorption
      res(ll,2) = tmp(0) - tmp(1); // cd ext
      res(ll,3) = tmp(2) - tmp(3); // cd abs
    }
    if(progress == 1)
      cout<<"\n";

    return res ;
  } 


arma::mat circular_dichroism_spectrum2(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
				      const arma::mat& Euler, const arma::mat& QMC, \
				      const int full, const int progress)
  {

    int N = kn.n_elem, Nr = R.n_rows, ll;
    // cout << N << "\n";
    arma::mat res(N,4);
    arma::cx_mat beta(3,Nr);
    arma::colvec tmp(4);
    arma::cx_mat A(3*Nr,3*Nr), polar(3*Nr,3*Nr);

    for(ll=0; ll<N; ll++){ // loop over kn   
      if(progress == 1)
	progress_bar(ll+1,N);
      beta = reshape(Beta.row(ll), 3, Nr, 1); 
      A = interaction_matrix(R, kn[ll], beta, Euler, full);
      polar = diagonal_polarisability(beta, Euler);
      tmp = averaging2(R, A, polar, kn[ll], QMC);

      res(ll,0) = 0.5*(tmp(0) + tmp(1)); // extinction 
      res(ll,1) = 0.5*(tmp(2) + tmp(3)); // absorption
      res(ll,2) = tmp(0) - tmp(1); // cd ext
      res(ll,3) = tmp(2) - tmp(3); // cd abs
    }
    if(progress == 1)
      cout<<"\n";

    return res ;
  } 


RCPP_MODULE(cd){
       using namespace Rcpp ;

       function( "circular_dichroism_spectrum2", &circular_dichroism_spectrum2, \
		 "Calculates the orientation-averaged CD spectrum for absorption and extinction" ) ;

       function( "circular_dichroism_spectrum", &circular_dichroism_spectrum, \
		 "Calculates the orientation-averaged CD spectrum for absorption and extinction" ) ;

}






