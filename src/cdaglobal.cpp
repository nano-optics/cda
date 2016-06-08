//
// These functions provide an alternative to cda.cpp
// whereby they update one argument, instead of returning a copy
// This is somewhat faster, but less transparent to use (hence the two versions)
//

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"
#include "cdaglobal.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Block-matrix of polarisabilities
//
// Alpha:  3N vector of principal polarizabilities
// Angles:  3xN matrix of cluster angles to apply
// updates a 3x3xN complex cube of polarizabilities
// [[Rcpp::export]]
void cpp_alpha_blocks_update(const arma::cx_colvec& Alpha,
                      const arma::mat& Angles,
                      arma::cx_cube& AlphaBlocks) {

  const int N = Angles.n_cols;
  arma::mat Rot(3,3);
  int ii=0;

  // loop over N particles
  for(ii=0; ii<N; ii++){

    Rot = cpp_euler_passive(Angles(0,ii), Angles(1,ii), Angles(2,ii));
    AlphaBlocks.slice(ii) =  Rot.st() * diagmat(Alpha.subvec(ii*3, ii*3+2)) * Rot;

  }

}


//
// Polarization from local fields and polarisabilities
//
// E: 3NxNi complex matrix of local field
// AlphaBlocks:  3x3xN block-matrix of polarizabilities
// updates the 3NxNi complex matrix of polarization
// [[Rcpp::export]]
void cpp_polarization_update(const arma::cx_mat& E,
                              const arma::cx_cube& AlphaBlocks,
                              arma::cx_mat& P) {

  const int N = AlphaBlocks.n_slices, Ni = E.n_cols;
  arma::cx_mat alphatmp(3,3);
  int ii=0;

  // loop over N particles
  for(ii=0; ii<N; ii++){
    alphatmp =  AlphaBlocks.slice(ii);
    P.submat(ii*3, 0, ii*3+2, Ni-1) = alphatmp * E.submat(ii*3, 0, ii*3+2, Ni-1);
  }

}

// Update the full interaction matrix
//
// R:  Nx3 matrix of positions
// kn:  incident wavenumber (scalar)
// AlphaBlocks:  3x3xN blocks of polarizabilities
// updates the 3Nx3N interaction matrix
// [[Rcpp::export]]
void cpp_interaction_matrix_update(const arma::mat& R, const double kn,
                                   const arma::cx_cube& AlphaBlocks,
                                   arma::cx_mat& A) {

  const int N = R.n_cols;

  // constants
  const arma::cx_double i = arma::cx_double(0,1);
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );

  // temporary vars

  int jj=0, kk=0;
  arma::mat rk_to_rj = arma::mat(3,1), rjkhat = arma::mat(3,1) ,
    rjkrjk = arma::mat(3,3);

  double rjk;
  arma::cx_mat Ajk = arma::cx_mat(3,3), alphajj = Ajk, alphakk=Ajk;

  // nested for loop over N dipoles
  for(jj=0; jj<N; jj++)
  {
    alphajj =  AlphaBlocks.slice(jj);

    for(kk=jj+1; kk<N; kk++)
    {
      alphakk =  AlphaBlocks.slice(kk);

      rk_to_rj = R.col(jj) - R.col(kk) ;
      rjk = norm(rk_to_rj, 2);
      rjkhat = rk_to_rj / rjk;
      rjkrjk =  rjkhat * rjkhat.st();

      Ajk = exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) +
        (i*kn*rjk - arma::cx_double(1,0)) / (rjk*rjk) * (3*rjkrjk - I3)) ;

      // assign block
      A.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Ajk * alphakk;
      // symmetric block
      A.submat(kk*3,jj*3,kk*3+2,jj*3+2) = Ajk.st() * alphajj;

    } // end kk
  } // end jj

}


// Update the propagator matrix
//
// R:  Nx3 matrix of positions
// kn:  incident wavenumber (scalar)
// AlphaBlocks:  3x3xN blocks of polarizabilities
// updates the 3Nx3N propagator matrix G
// [[Rcpp::export]]
void cpp_propagator_update(const arma::mat& R, const double kn,
                                   const arma::cx_cube& AlphaBlocks,
                                   arma::cx_mat& G) {
  
  const int N = R.n_cols;
  
  // constants
  const arma::cx_double i = arma::cx_double(0,1);
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );
  
  // temporary vars
  
  int jj=0, kk=0;
  arma::mat rk_to_rj = arma::mat(3,1), rjkhat = arma::mat(3,1) ,
    rjkrjk = arma::mat(3,3);
  
  double rjk;
  arma::cx_mat Gjk = arma::cx_mat(3,3), alphajj = Gjk, alphakk=Gjk;
  
  // nested for loop over N dipoles
  for(jj=0; jj<N; jj++)
  {
    alphajj =  AlphaBlocks.slice(jj);
    
    for(kk=jj+1; kk<N; kk++)
    {
      alphakk =  AlphaBlocks.slice(kk);
      
      rk_to_rj = R.col(jj) - R.col(kk) ;
      rjk = norm(rk_to_rj, 2);
      rjkhat = rk_to_rj / rjk;
      rjkrjk =  rjkhat * rjkhat.st();
      
      Gjk = -exp(i*kn*rjk) / rjk *  (kn*kn*(rjkrjk - I3) +
        (i*kn*rjk - arma::cx_double(1,0)) / (rjk*rjk) * (3*rjkrjk - I3)) ;
      
      // assign block
      G.submat(jj*3,kk*3,jj*3+2,kk*3+2) = Gjk * alphakk;
      // symmetric block
      G.submat(kk*3,jj*3,kk*3+2,jj*3+2) = Gjk.st() * alphajj;
      
    } // end kk
  } // end jj
  
}
