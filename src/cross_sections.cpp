// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>
#include "cross_sections.h"
#include "utils.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Calculate the extinction cross-section for multiple incident angles
//
// kn: scalar wavenumber
// P: complex 3NxNi matrix of polarization
// Ein: complex 3NxNs matrix of incident fields
// return the extinction cross-section for all incident angles
// [[Rcpp::export]]
arma::colvec cpp_extinction(const double kn,
                            const arma::cx_mat& P,
                            const arma::cx_mat& Ein)
{
  int Ni = P.n_cols, ii=0;
  arma::colvec cext(Ni);

  for (ii=0; ii<Ni; ii++)
  {
    cext(ii) = imag(cdot(Ein.col(ii), P.col(ii)));
  }

  return  4*arma::datum::pi*kn*cext;
}


//
// Calculate the absorption cross-section for multiple incident angles
//
// kn: scalar wavenumber
// P: complex 3NxNi matrix of polarization
// E: complex 3NxNi matrix of local fields
// return the absorption cross-section for all incident angles
// [[Rcpp::export]]
arma::colvec cpp_absorption(const double kn,
                            const arma::cx_mat& P,
                            const arma::cx_mat& E)
{
  int Ni = P.n_cols, ii=0;
  arma::colvec cabs(Ni);

  for (ii=0; ii<Ni; ii++)
  {
    cabs(ii) = imag(cdot(E.col(ii), P.col(ii))) -
      kn*kn*kn* 2/3 * real(cdot(P.col(ii), P.col(ii)));
  }
  return  4*arma::datum::pi*kn*cabs;

}

//
// Calculate the scattering cross-section for multiple incident angles
//
// R: 3xN matrix of dipole positions
// angles: 2xNs matrix of scattering angles
// weights: Ns vector of scattering quadrature weights
// kn: scalar wavenumber
// P: complex 3NxNi matrix of polarization
// return the scattering cross-section for all incident angles
// [[Rcpp::export]]
arma::colvec cpp_scattering(const arma::mat& R,
                            const arma::mat& angles,
                            const arma::colvec& weights,
                            const double kn,
                            const arma::cx_mat& P)
{
  const int Nr = R.n_cols, Ni = P.n_cols, Ns = angles.n_cols;
  int ii=0, jj=0;

  arma::cx_mat G = arma::zeros<arma::cx_mat>(3,3);
	const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );

  arma::mat Rm(3,3);

  arma::cx_mat Esca(3, Ni);
  arma::mat nn(3,3);
  arma::colvec n(3), rj(3);
  double nrj;
	const arma::cx_double i = arma::cx_double(0,1);
	arma::cx_double phase;

  arma::mat Isca(Ns, Ni);
  arma::rowvec csca(Ni);


  for (ii=0; ii<Ns; ii++) // loop over scattering angles
  {

    Rm = cpp_euler_active(angles(0,ii), angles(1,ii), angles(2,ii));
    n = Rm.col(2);
    nn = n*n.st();
    G = (I3 - nn);
    Esca = arma::zeros<arma::cx_mat>(3, Ni);
    for (jj=0; jj<Nr; jj++) // loop over dipoles
    {
      rj = R.col(jj);
      nrj = dot(n, rj);
      phase = exp(-i*kn*nrj) ;
      Esca = Esca + phase * G * P.submat(jj*3, 0, jj*3+2, Ni-1);
    }

    //  Rcpp::Rcout << "here" << "\n";
    // Esca is now the net FF in direction ii
    Isca.row(ii) = real(sum(Esca % conj(Esca), 0)); // |Esca|^2

    //Rcpp::Rcout << Isca.row(ii) << "\n";

  }

  //Rcpp::Rcout << "then" << "\n";
  csca =  weights.st() * Isca;
  // Rcpp::Rcout << csca << "\n";

  return  pow(kn,4) * 4*arma::datum::pi*csca.st();

}
