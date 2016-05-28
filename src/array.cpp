// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "array.h"

//
// Array factor for infinite square array at normal incidence
//
// kn: vector of incident wavenumbers
// R: Nx3 matrix of positions
// [[Rcpp::export]]
arma::cx_colvec cpp_array_factor(const arma::colvec& kn, const arma::mat& R) {

  const int Nk = kn.n_elem, Nr = R.n_rows;
  const arma::cx_double i = arma::cx_double(0,1);

  double kk,rmod,kr,sin2theta;
  arma::cx_double A,B,S;
  arma::cx_colvec array(Nk);

  int ik=0, iR=0; // loop over wavelengths
  for(ik=0; ik<Nk; ik++){
    kk = kn(ik);
    S = 0;
    //loop over R rows
    for(iR = 0; iR < Nr; iR++){
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
