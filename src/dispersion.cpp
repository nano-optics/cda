// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"
#include "cdaglobal.h"
#include "dispersion.h"
#include "cg.h"
#include "iterative.h"
#include "incident.h"
#include "cross_sections.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Angle-resolved cross-sections for multiple directions of incidence
//
// R: Nx3 matrix of positions
// A: 3Nx3N interaction matrix
// AlphaBlocks:  3x3xN blocks of polarizabilities
// Alpha: 3NxNl matrix of principal polarisabilities
// kn:  incident wavenumber (scalar)
// medium:  incident medium (scalar)
// Incidence:  Ni vector of incident beam directions
// Axes:  Ni vector of incident beam axes
// ScatteringNodes:  2xNs matrix of scattered field directions
// ScatteringWeights:  Ns vector of scattered field quadrature weights
// polarisation: integer flag to switch between linear and circular polarisation
// inversion: integer flag to use different inversion methods
// maxiter:  max number of iterations
// tol:  tolerance
// [[Rcpp::export]]
arma::mat cpp_dispersion(const arma::mat& R,
  const arma::cx_mat& A,
  const arma::cx_cube& AlphaBlocks,
  const double kn, const double medium,
  const arma::vec& Incidence,
  const arma::ivec& Axes,
  const arma::mat& ScatteringNodes,
  const arma::colvec& ScatteringWeights,
  const int polarisation,
  const int inversion,
  const int maxiter, const double tol)
  {
    const int N = R.n_cols, Ni = Incidence.n_elem;

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

    Ein.submat(0,0,3*N-1,Ni-1) = cpp_incident_field_axis(Evec1, kn*khat, R, Incidence, Axes);
    Ein.submat(0,Ni,3*N-1,2*Ni-1) = cpp_incident_field_axis(Evec2, kn*khat, R, Incidence, Axes);

    arma::mat res(Ni, 6) ;
    arma::colvec xsec(2*Ni); // temporary storage of cross-sections

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

      // return angle-dependent cext, cabs, csca
      // cext
      xsec = cpp_extinction(kn, P, Ein);
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

      //res.col(2) = res.col(0) - res.col(1);// first polar
      //res.col(5) = res.col(3) - res.col(4);// second polar

      return res ;
    }



    //
    // Angle-resolved spectra for linear or circular polarisations
    //
    // kn: vector of incident wavenumbers
    // medium: surrounding refractive index
    // R: Nx3 matrix of positions
    // Alpha: 3NxNl matrix of principal polarisabilities
    // Angles:  Nx3 matrix of particle rotation angles
    // Incidence:  Ni vector of incident beam directions
    // Axes:  Ni vector of incident beam axes
    // ScatteringNodes:  2xNs matrix of scattered field directions
    // ScatteringWeights:  Ns vector of scattered field quadrature weights
    // polarisation: integer flag to switch between linear and circular polarisation
    // inversion: integer flag to use different inversion methods
    // maxiter:  max number of iterations
    // tol:  tolerance
    // progress: logical flag to display progress bars
    // [[Rcpp::export]]
    arma::cube cpp_dispersion_spectrum(const arma::colvec kn,
      const double medium,
      const arma::mat& R,
      const arma::cx_mat& Alpha,
      const arma::mat& Angles,
      const arma::vec& Incidence,
      const arma::ivec& Axes,
      const arma::mat& ScatteringNodes,
      const arma::vec& ScatteringWeights,
      const int polarisation,
      const int inversion,
      const int maxiter,
      const double tol,
      const bool progress)
      {

        const int Ni = Incidence.n_elem;
        int Nl = kn.n_elem, Nr = R.n_cols, ll;
        arma::cube res(Ni, 6, Nl);
        arma::mat tmp(Ni, 6);

        // global, will update at each wavelength
        arma::cx_cube AlphaBlocks = arma::zeros<arma::cx_cube>(3, 3, Nr);

        // global, will update at each wavelength
        arma::cx_mat A(3*Nr, 3*Nr );
        if(inversion < 2) { // full interaction matrix
          A.eye();
        } else if (inversion == 2){ // only propagator G = I - A
          A.zeros();
        }

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

          tmp = cpp_dispersion(R, A, AlphaBlocks, kn(ll), medium,
          Incidence, Axes, ScatteringNodes, ScatteringWeights, polarisation, inversion, maxiter, tol);
          // Rcpp::Rcout << tmp << "\n";

          res.slice(ll) = 1.0/Nr * tmp;
        }
        if(progress)
        Rcpp::Rcout << "\n";

        return res ;
      }
