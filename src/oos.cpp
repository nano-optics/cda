// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"
#include "cda.h"
#include "oos.h"
#include "incident.h"
#include "cross_sections.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;


/////////////////////////////////
// order of scattering variant //
/////////////////////////////////


// iterate the local field and polarisation
//
// R is the Nx3 matrix of positions
// kn: scalar wavenumber
// Ein is the 3NxNi incident field
// AlphaBlocks is the 3x3xN block-matrix of polarizabilities
// E is the 3NxNi local field
// P: complex 3NxNi matrix of polarization
// return current extinction cross-sections (for convergence check)
// side-effect: update P and E
// [[Rcpp::export]]
arma::colvec cpp_iterate_field(const arma::mat& R,
                               const double kn, const double medium,
                               const arma::cx_mat& Ein,
                               const arma::cx_cube& AlphaBlocks,
                               arma::cx_mat& Eloc,
                               arma::cx_mat& P){

  const int N = R.n_cols;
  const int Ni = Ein.n_cols;

  // tmp variables
  const arma::cx_mat I3 = arma::eye<arma::cx_mat>( 3, 3 );
  const arma::cx_double i = arma::cx_double(0,1);
  arma::cx_mat Gjk = arma::cx_mat(3,3);
  double rjk;
  arma::mat rk_to_rj = arma::mat(3,1);
  arma::mat rjkhat   = arma::mat(3,1);
  arma::mat rjkrjk   = arma::mat(3,3);
  int ll=0, jj=0, kk=0;
  arma::colvec cext(Ni); // extinction for convergence test
  arma::cx_vec Edip(3);
  Eloc = Ein; // incident field +
  // all pairwise dipole fields
  for(jj=0; jj<N; jj++)
  {
    for(kk=jj+1; kk<N; kk++)
    {
      rk_to_rj = R.col(jj) - R.col(kk) ;
      rjk = norm(rk_to_rj, 2);
      rjkhat = rk_to_rj / rjk;
      rjkrjk =  rjkhat * rjkhat.st();
      // 3x3 propagator
      Gjk = exp(i*kn*rjk) / rjk *  (kn*kn*(I3 - rjkrjk) +
        (i*kn*rjk - arma::cx_double(1,0)) /
                                      (rjk*rjk) * (I3 - 3*rjkrjk)) ;

      // update E = Ein + GP
      // where P is from a previous iteration
      // loop over IncidenceNodes angles (columns of E)
      for(ll=0; ll<Ni; ll++){
        // field of dipole kk evaluated at jj
        Edip = Gjk * P.submat(kk*3, ll, kk*3+2, ll);
        Eloc.submat(jj*3, ll, jj*3+2, ll) += Edip;
        // field of dipole jj evaluated at kk
        Edip = Gjk.st() * P.submat(jj*3, ll, jj*3+2, ll);
        Eloc.submat(kk*3, ll, kk*3+2, ll) += Edip;
      }
    }
  }
  // update the polarization
  P = cpp_polarization(Eloc, AlphaBlocks);
  cext = cpp_extinction(kn, P, Ein);

  return(cext);
}


// Performs full angular averaging for both circular polarisations
//
// R is the Nx3 matrix of positions
// AlphaBlocks is the 3Nx3 block-matrix of polarizabilities
// kn is the incident wavenumber (scalar)
// Angles is the Nix3 matrix of incident beam angles
// Weigths is the Ni vector of quadrature IncidenceWeights
// maxiter is the max number of  iterations
// tol is the tolerance
// [[Rcpp::export]]
arma::rowvec cpp_oa_oos(const arma::mat& R,
                       const arma::cx_cube& AlphaBlocks,
                       const double kn, const double medium,
                       const arma::mat& IncidenceNodes,
                       const arma::colvec& IncidenceWeights,
                       const arma::mat& ScatteringNodes,
                       const arma::colvec& ScatteringWeights,
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
  arma::cx_mat Ein(3*N, 2*Ni), Eloc=Ein, P=Ein;
  const arma::colvec  khat="0 0 1;";

  arma::cx_colvec Evec = "(0,1) (1,0) (0,0);";
  Evec = arma::datum::sqrt2/2 * Evec ;
  Ein.submat(0,0,3*N-1,Ni-1) = cpp_incident_field(Evec, kn*khat, R, IncidenceNodes);
  Evec = "(1,0) (0,1) (0,0);";
  Evec = arma::datum::sqrt2/2 * Evec ;
  Ein.submat(0,Ni,3*N-1,2*Ni-1) = cpp_incident_field(Evec, kn*khat, R, IncidenceNodes);

  // tmp variables
  int iter=0;
  double rel_error=1e10;
  arma::colvec xsec = zeros<arma::colvec>(2*Ni), tmp=xsec;

  // -----------------------
  // set up
  // -----------------------
  iter=0;
  rel_error=1e10;
  Eloc=Ein;
  P = cpp_polarization(Eloc, AlphaBlocks);
  tmp = cpp_extinction(kn, P, Eloc);

  // order-of-scattering iterations
  while((iter < maxiter) && (rel_error > tol)){
    xsec = cpp_iterate_field(R, kn, medium, Ein, AlphaBlocks, Eloc, P);
    //Rcpp::Rcout << xsec << "\n";
    // Note E and P have been updated
    rel_error = max(abs((xsec - tmp) / (xsec + tmp)));
    tmp = xsec;
    iter++;
  }
  // xsec contains the latest iteration of extinction
  res(0) = 0.5*dot(xsec, IncidenceWeights1); // cext
  res(3) = dot(xsec, IncidenceWeights2);     // dext
  xsec = cpp_absorption(kn, P, Eloc);
  res(1) = 0.5*dot(xsec, IncidenceWeights1); // cabs
  res(4) = dot(xsec, IncidenceWeights2);     // dabs

  // csca
  xsec =  cpp_scattering(R, ScatteringNodes, ScatteringWeights, kn, P);
  res(2) = 0.5*dot(xsec, IncidenceWeights1); // csca
  res(5) = dot(xsec, IncidenceWeights2);    //  dsca

  return res ;
}


// Angle-averaged spectra for LCP and RCP polarisations
//
// kn is the vector of incident wavenumbers
// R is the Nx3 matrix of positions
// Alpha is the 3N vector of principal polarisabilities
// Euler is the Nx3 matrix of particle rotation angles
// Angles is the Nix3 matrix of incident beam angles
// Weigths is the Ni vector of quadrature IncidenceWeights
// full is a logical flag to switch off retardation terms
// maxiter is the max number of cg iterations
// tol is the cg tolerance
// progress is a logical flag to display progress bars
// [[Rcpp::export]]
arma::mat cpp_oa_spectrum_oos(const arma::colvec kn,
                             const double medium,
                             const arma::mat& R,
                             const arma::cx_mat& Alpha,
                             const arma::mat& Angles,
                             const arma::mat& IncidenceNodes,
                             const arma::colvec& IncidenceWeights,
                             const arma::mat& ScatteringNodes,
                             const arma::colvec& ScatteringWeights,
                             const bool full,
                             const int maxiter,
                             const double tol,
                             const bool progress)
{

  int Nl = kn.n_elem, Nr = R.n_cols, ll;

  arma::mat res(Nl,6);
  arma::rowvec tmp(6);
  arma::cx_cube AlphaBlocks(3,3,Nr);

  for(ll=0; ll<Nl; ll++){ // loop over kn
    if(progress)
      progress_bar(ll+1,Nl);

    AlphaBlocks = cpp_alpha_blocks(Alpha.col(ll), Angles);
    tmp = cpp_oa_oos(R, AlphaBlocks, kn(ll), medium,
                     IncidenceNodes, IncidenceWeights,
                     ScatteringNodes, ScatteringWeights, maxiter, tol);
    res.row(ll) = 1.0/Nr * tmp;
  }
  if(progress)
    Rcpp::Rcout << "\n";

  return res ;
}

//////////////////////////////////////////////////

arma::mat cpp_dispersion_oos(const arma::mat& R,
                            const arma::cx_cube& AlphaBlocks,
                            const double kn, const double medium,
                            const arma::vec& IncidenceNodes,
                            const arma::ivec& IncidenceAxes,
                            const arma::mat& ScatteringNodes,
                            const arma::colvec& ScatteringWeights,
                            const int polarisation,
                            const int maxiter, const double tol)
{
    const int N = R.n_cols, Ni = IncidenceNodes.n_elem;
    arma::mat res(Ni, 6) ;

    // incident field
    arma::cx_mat Ein(3*N, 2*Ni), Eloc(3*N, 2*Ni);
    arma::cx_mat P(3*N, 2*Ni);
    const arma::colvec  khat="0 0 1;";

  arma::cx_colvec Evec1, Evec2;
  if(polarisation == 0){ // linear
    Evec1="(1,0) (0,0) (0,0);", Evec2="(0,0) (1,0) (0,0);";
  } else { // circular
    Evec1="(0,1) (1,0) (0,0);", Evec2="(1,0) (0,1) (0,0);";
    Evec1 = arma::datum::sqrt2/2 * Evec1 ;
    Evec2 = arma::datum::sqrt2/2 * Evec2 ;
  }

  Ein.submat(0,0,3*N-1,Ni-1) = cpp_incident_field_axis(Evec1, kn*khat, R, IncidenceNodes, IncidenceAxes);
  Ein.submat(0,Ni,3*N-1,2*Ni-1) = cpp_incident_field_axis(Evec2, kn*khat, R, IncidenceNodes, IncidenceAxes);

  // tmp variables
  int iter=0;
  double rel_error=1e10;
  arma::colvec xsec = zeros<arma::colvec>(2*Ni), tmp=xsec;
  // -----------------------
  // set up
  // -----------------------
  iter=0;
  rel_error=1e10;
  Eloc=Ein;
  P = cpp_polarization(Eloc, AlphaBlocks);
  tmp = cpp_extinction(kn, P, Eloc);

  // order-of-scattering iterations
  while((iter < maxiter) && (rel_error > tol)){
    xsec = cpp_iterate_field(R, kn, medium, Ein, AlphaBlocks, Eloc, P);
    //Rcpp::Rcout << xsec << "\n";
    // Note E and P have been updated
    rel_error = max(abs((xsec - tmp) / (xsec + tmp)));
    tmp = xsec;
    iter++;
  }

  // return angle-dependent cext, cabs, csca
  // cext
  // xsec contains the latest iteration of extinction
  res.col(0) = xsec.subvec(0, Ni-1);    // first polar
  res.col(3) = xsec.subvec(Ni, 2*Ni-1); // second polar
  // cabs
  xsec = cpp_absorption(kn, P, Eloc);
  res.col(1) = xsec.subvec(0, Ni-1);    // first polar
  res.col(4) = xsec.subvec(Ni, 2*Ni-1); // second polar

  // csca
  xsec =  cpp_scattering(R, ScatteringNodes, ScatteringWeights, kn, P);
  res.col(2) = xsec.subvec(0, Ni-1);    // first polar
  res.col(5) = xsec.subvec(Ni, 2*Ni-1); // second polar


  return res ;
}


//
// Angle-resolved spectra for linear or circular polarisations
//
// kn is the vector of incident wavenumbers
// Alpha is the 3N vector of principal polarisabilities
// R is the Nx3 matrix of positions
// Angles is the Nx3 matrix of particle rotation IncidenceNodes
// IncidenceNodes is the Ni vector of incident beam IncidenceNodes
// IncidenceAxes is the Ni vector of incident beam IncidenceAxes
// polarisation is an integer flag to switch between linear and circular polarisation
// full is a logical flag to switch off retardation terms
// cg logical flag to use conjugate gradient solver
// born logical flag, use first Born approx for cg solver
// maxiter is the max number of cg iterations
// tol is the cg tolerance
// progress is a logical flag to display progress bars
// [[Rcpp::export]]
arma::cube cpp_dispersion_spectrum_oos(const arma::colvec kn,
             const double medium,
			       const arma::mat& R,
			       const arma::cx_mat& Alpha,
			       const arma::mat& Angles,
			       const arma::vec& IncidenceNodes,
			       const arma::ivec& IncidenceAxes,
             const arma::mat& ScatteringNodes,
             const arma::colvec& ScatteringWeights,
			       const int polarisation,
             const bool full,
             const int maxiter,
             const double tol,
             const bool progress)
  {

    const int Ni = IncidenceNodes.n_elem;
    int Nl = kn.n_elem, Nr = R.n_cols, ll;
    arma::cube res(Ni, 6, Nl);
    arma::mat tmp(Ni, 6);
    arma::cx_cube AlphaBlocks(3,3,Nr);

    for(ll=0; ll<Nl; ll++){ // loop over kn
      if(progress)
	progress_bar(ll+1,Nl);

      AlphaBlocks = cpp_alpha_blocks(Alpha.col(ll), Angles);
      tmp = cpp_dispersion_oos(R, AlphaBlocks, kn(ll), medium,
      IncidenceNodes, IncidenceAxes, ScatteringNodes, ScatteringWeights, polarisation, maxiter, tol);
      // Rcpp::Rcout << tmp << "\n";

      res.slice(ll) = 1.0/Nr * tmp;
    }
    if(progress)
      Rcpp::Rcout << "\n";

    return res ;
  }
