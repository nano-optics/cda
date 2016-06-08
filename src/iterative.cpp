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
int cpp_iterate_field(const arma::cx_mat& E0,
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
  while((error > tol) && (niter < maxiter)){
    
    Etmp = G * Etmp; // new order of scattering
    Eloc = Eloc + Etmp; // update total field
    cpp_polarization_update(Eloc, AlphaBlocks, P);// update P
    
    cextnew = cpp_extinction(kn, P, E0); // testing convergence
    error = max(abs((cextnew - cext)) / abs(cextnew + cext));
    cext = cextnew;
    niter++;
  }
  
  return(niter);
}

