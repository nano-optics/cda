// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
#include "cda.h"
#include "iterative.h"
#include "cdaglobal.h"
#include "incident.h"
#include "cross_sections.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;


// iterate the local field and polarisation
//
// R is the Nx3 matrix of positions
// kn: scalar wavenumber
// Ein is the 3NxNi incident field
// AlphaBlocks is the 3x3xN block-matrix of polarizabilities
// E is the 3NxNi local field
// P: complex 3NxNi matrix of polarization
// side-effect: update E and P
bool cpp_iterate_field(const arma::cx_mat& E0,
                      const arma::cx_mat& G,
                      const arma::cx_cube& AlphaBlocks,
                      const double kn,
                      const double tol,
                      const int maxiter,
                      arma::cx_mat& Eloc,
                      arma::cx_mat& P){

  // initialise
  arma::cx_mat Etmp = E0; // First Born Approximation
  Eloc = Etmp;
  cpp_polarization_update(Eloc, AlphaBlocks, P); // update P
  arma::colvec cext, cextnew;
  cext = cpp_extinction(kn, P, E0); // testing convergence

  double error = 1e6; // init value doesn't matter as long as > tol
  int niter = 1;
  bool expectation = true; // still work to do to reach target
  while(expectation && (niter < maxiter)){

    Etmp = G * Etmp; // new order of scattering
    Eloc = Eloc + Etmp; // update total field
    cpp_polarization_update(Eloc, AlphaBlocks, P);// update P

    cextnew = cpp_extinction(kn, P, E0); // testing convergence
    error = max(abs((cextnew - cext)) / abs(cextnew + cext));
    cext = cextnew;
    niter++;
    // still work to do to reach target?
    expectation = (error > tol);
  }
  if(expectation){
    Rcpp::warning("O-O-S target not reached: rel. err. %.2e in %i iter", error, niter);
  }
  return(expectation);
}
