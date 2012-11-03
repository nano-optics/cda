#include "array.h"

#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;
using namespace arma;

cx_colvec array_factor(const colvec& k, const mat& R) {
  
  const int Nk = k.n_elem, NR = R.n_rows;
  const cx_double i = cx_double(0,1);

  double kk,rmod,kr,sin2theta;
  cx_double A,B,S;
  cx_colvec array(Nk);

  int ik=0, iR=0; // loop over wavelengths
  for(ik=0; ik<Nk; ik++){
    kk = k(ik);
    S = 0;
    //loop over R rows
    for(iR = 0; iR < NR; iR++){
      rmod = norm(R.row(iR), 2);
      kr = kk*rmod;
      sin2theta = pow(R(iR,1) / rmod, 2);
      A = (2 - 3*sin2theta) / pow(kr,3);    
      B = -i*kr*A;
      S += pow(kk,3) * exp(i*kr) * (A + B + sin2theta/kr);
    } // done all R
    
    array(ik) = S;
    
  } // done all wavelengths
  

  return array;
}



RCPP_MODULE(array){
       using namespace Rcpp ;
       function( "array_factor", &array_factor ) ;
}


