// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>

#include "incident.h"
#include "utils.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

//
// Calculate the incident field at each dipole location
//
// for multiple euler angles of incidence
// E0 is the normalised incident electric field before rotation
// k is the incident wavevector before rotation
// R is the 3xN matrix of positions
// Angles is the 3xNi matrix of incident beam angles
// returns a 3NxNi complex matrix of incident fields
// where each column corresponds to a specific incident angle
// [[Rcpp::export]]
arma::cx_mat cpp_incident_field(const arma::cx_colvec& Evec,
			    const arma::colvec& kvec,
			    const arma::mat& R,
			    const arma::mat& Incidence)
{
  const int Ni = Incidence.n_cols;
  const int N = R.n_cols;
  const arma::cx_double i = arma::cx_double(0,1);
  arma::mat Rot(3,3);
  arma::cx_mat Ein = arma::cx_mat(3*N,Ni);
  arma::cx_colvec Evec_r(3);
  arma::colvec kvec_r(3);
  arma::rowvec kR(N);
  arma::cx_rowvec expikr(N);
  int jj=0;
  for(jj=0; jj<Ni; jj++)
    {
      Rot = cpp_euler_active(Incidence(0,jj), Incidence(1,jj), Incidence(2,jj));
      kvec_r = Rot * kvec;
      Evec_r = Rot * Evec;
      kR = strans(kvec_r) * R; // [1 x N]
      expikr = exp(i * kR); // [1 x N]

      Ein.col(jj) = vectorise(Evec_r * expikr, 0);
    }
  return(Ein);
}

//
// Incident field along multiple axes
//
// E0 is the normalised electric field
// k is the wavevector
// R is the Nx3 matrix of positions
// Axes is a vector of integer codes corresponding to x, y, z
// Angles is a vector of rotation angles around Axes
// [[Rcpp::export]]
arma::cx_mat cpp_incident_field_axis(const arma::cx_colvec& Evec,
			     const arma::colvec& kvec,
			     const arma::mat& R,
			     const arma::colvec& Incidence,
					 const arma::ivec& Axes)
{
  const int Ni = Incidence.n_elem;
  const int N = R.n_cols;
  const arma::cx_double i = arma::cx_double(0,1);
  arma::mat Rot(3,3);
  arma::cx_mat Ein = arma::cx_mat(3*N,Ni);
  arma::cx_colvec Evec_r(3);
	arma::colvec kvec_r(3);
  arma::rowvec kR(N);
  arma::cx_rowvec expikr(N);
  int jj=0;
  for(jj=0; jj<Ni; jj++)
    {
      Rot = cpp_axis_rotation(Incidence(jj), Axes(jj));
			kvec_r = Rot * kvec;
      Evec_r = Rot * Evec;
      kR = strans(kvec_r) * R; // [1 x N]
      expikr = exp(i * kR); // [1 x N]

      Ein.col(jj) = vectorise(Evec_r * expikr, 0);

    }

  return(Ein);
}
