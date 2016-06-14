// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"
#include "cdaglobal.h"
#include "averaging.h"
#include "cg.h"
#include "iterative.h"
#include "incident.h"
#include "cross_sections.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;

//
// Performs full angular averaging for both circular polarisations
//
// R: Nx3 matrix of positions
// A: 3Nx3N interaction matrix
// AlphaBlocks:  3x3xN blocks of polarizabilities
// Alpha: 3NxNl matrix of principal polarisabilities
// kn:  incident wavenumber (scalar)
// medium:  incident medium (scalar)
// IncidenceNodes:  2xNi matrix of incident field directions
// IncidenceWeights:  Ni vector of incident field quadrature weights
// ScatteringNodes:  2xNs matrix of scattered field directions
// ScatteringWeights:  Ns vector of scattered field quadrature weights
// inversion: integer flag to use different inversion methods
// maxiter:  max number of iterations
// tol:  tolerance
// [[Rcpp::export]]
arma::rowvec cpp_oa(const arma::mat& R,
                    const arma::cx_mat& A,
                    const arma::cx_cube& AlphaBlocks,
                    const double kn,
                    const double medium,
                    const arma::mat& IncidenceNodes,
                    const arma::colvec& IncidenceWeights,
                    const arma::mat& ScatteringNodes,
                    const arma::colvec& ScatteringWeights,
                    const int inversion,
                    const int maxiter, const double tol)
{

  const int N = R.n_cols, Ni = IncidenceNodes.n_cols;
  // combine IncidenceWeights for both polarisations
  // to get average and dichroism
  arma::colvec IncidenceWeights1(2*Ni), IncidenceWeights2(2*Ni);
  IncidenceWeights1 = join_vert(IncidenceWeights,IncidenceWeights);
  IncidenceWeights2 = join_vert(IncidenceWeights,-1*IncidenceWeights);
  arma::rowvec res(6) ;

  // incident field
  arma::cx_mat Ein(3*N, 2*Ni), Eloc(3*N, 2*Ni), P(3*N, 2*Ni);
  const arma::colvec  khat="0 0 1;";

  arma::cx_colvec Evec = "(0,1) (1,0) (0,0);";
  Evec = arma::datum::sqrt2/2 * Evec ;
  Ein.submat(0,0,3*N-1,Ni-1) = cpp_incident_field(Evec, kn*khat, R, IncidenceNodes);
  Evec = "(1,0) (0,1) (0,0);";
  Evec = arma::datum::sqrt2/2 * Evec ;
  Ein.submat(0,Ni,3*N-1,2*Ni-1) = cpp_incident_field(Evec, kn*khat, R, IncidenceNodes);

  arma::colvec xsec(Ni); // temporary storage of cross-sections

  if(inversion == 0) { // standard method
    Eloc = solve(A, Ein);
    cpp_polarization_update(Eloc, AlphaBlocks, P);
  } else if (inversion == 1){ // CG inversion
    arma::cx_mat guess = arma::zeros<arma::cx_mat>(3*N,2*Ni);
    Eloc = cpp_cg_solve(A, Ein, guess, maxiter, tol);
    cpp_polarization_update(Eloc, AlphaBlocks, P);
  } else if (inversion == 2){ // OOS solution (no inversion)
    bool expectation;
    // note: A is actually G here, not the same as above cases
    // we solve (I-G)E=Einc with G, G^2, G^3 etc.
    expectation = cpp_iterate_field(Ein, A, AlphaBlocks, kn,
                                    tol, maxiter, Eloc, P);
    // Eloc and P have now been updated by OOS
  }


  // cross section for all quadrature angles
  // averaged for polarisation and dichroism
  // cext cabs csca dext dabs dsca
  xsec = cpp_extinction(kn, P, Ein);

  res(0) = 0.5*dot(xsec, IncidenceWeights1); // cext
  res(3) = dot(xsec, IncidenceWeights2);     // dext
  xsec = cpp_absorption(kn, P, Eloc);
  res(1) = 0.5*dot(xsec, IncidenceWeights1); // cabs
  res(4) = dot(xsec, IncidenceWeights2);     // dabs

  // csca
  xsec =  cpp_scattering(R, ScatteringNodes, ScatteringWeights, kn, P);
  res(2) = 0.5*dot(xsec, IncidenceWeights1); // csca
  res(5) = dot(xsec, IncidenceWeights2);    //  dsca

	//Rcpp::Rcout << res << "\n";
  return res ;
}

//
// Angle-averaged spectra for LCP and RCP polarisations
//
// kn: vector of incident wavenumbers
// medium: surrounding refractive index
// R: Nx3 matrix of positions
// Alpha: 3NxNl matrix of principal polarisabilities
// Angles:  Nx3 matrix of particle rotation angles
// IncidenceNodes:  2xNi matrix of incident field directions
// IncidenceWeights:  Ni vector of incident field quadrature weights
// ScatteringNodes:  2xNs matrix of scattered field directions
// ScatteringWeights:  Ns vector of scattered field quadrature weights
// inversion: integer flag to use different inversion methods
// maxiter:  max number of iterations
// tol:  tolerance
// progress: logical flag to display progress bars
// [[Rcpp::export]]
arma::mat cpp_oa_spectrum(const arma::colvec kn,
                          const double medium,
                          const arma::mat& R,
                          const arma::cx_mat& Alpha,
                          const arma::mat& Angles,
                          const arma::mat& IncidenceNodes,
                          const arma::colvec& IncidenceWeights,
                          const arma::mat& ScatteringNodes,
                          const arma::colvec& ScatteringWeights,
                          const int inversion,
                          const int maxiter,
                          const double tol,
                          const bool progress)
{

  int Nl = kn.n_elem, Nr = R.n_cols, ll;

  arma::mat res(Nl,6);
  arma::rowvec tmp(6);

  // global, will update at each wavelength
  arma::cx_mat A(3*Nr, 3*Nr );
  if(inversion < 2) { // full interaction matrix
    A.eye();
  } else if (inversion == 2){ // only propagator G = I - A
    A.zeros();
  }

  arma::cx_cube AlphaBlocks = arma::zeros<arma::cx_cube>(3, 3, Nr);

  for(ll=0; ll<Nl; ll++){ // loop over kn
    if(progress)
      progress_bar(ll+1,Nl);

    // update in place
    cpp_alpha_blocks_update(Alpha.col(ll), Angles, AlphaBlocks);

    // if solve or cg, need full A, otherwise just G
    if(inversion < 2) { // full interaction matrix
      // update in place
      cpp_interaction_matrix_update(R, kn(ll), AlphaBlocks, A);
    } else if (inversion == 2){ // only propagator G = I - A
      // update in place
      cpp_propagator_update(R, kn(ll), AlphaBlocks, A);
    }
    // Rcpp::Rcout << A << "\n";

    tmp = cpp_oa(R, A, AlphaBlocks, kn(ll), medium,
                IncidenceNodes, IncidenceWeights, ScatteringNodes, ScatteringWeights,
                inversion, maxiter, tol);
    res.row(ll) = 1.0/Nr * tmp;
  }
  if(progress)
    Rcpp::Rcout << "\n";

  return res ;
}
