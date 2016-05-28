// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <iostream>

#include "utils.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace std;

void progress_bar(double x, double N)
  {
    // how wide you want the progress meter to be
    int totaldotz=40;
    double fraction = x / N;
    // part of the progressmeter that's already "full"
    int dotz = round(fraction * totaldotz);

    // create the "meter"
    int ii=0;
     Rprintf("%3.0f%% [",fraction*100);
    // part  that's full already
    for ( ; ii < dotz;ii++) {
      Rcpp::Rcout << "=";
    }
    // remaining part (spaces)
    for ( ; ii < totaldotz;ii++) {
       Rprintf(" ");
    }
    // and back to line begin -
    // do not forget the fflush to avoid output buffering problems!
    Rprintf("]\r");
    // fflush(stdout);
  }

// Euler rotation matrix
// ZYZ convention - active rotation
// to rotate a colvec V, do R * V
// [[Rcpp::export]]
arma::mat cpp_euler_active(const double phi, const double theta, const double psi)
  {
    arma::mat Rot(3,3);
    const double cosphi = cos(phi), cospsi = cos(psi), costheta = cos(theta);
    const double sinphi = sin(phi), sinpsi = sin(psi), sintheta = sin(theta);

    // note: indices have been simply swapped from passive version
    Rot(0,0) = cosphi*costheta*cospsi - sinphi*sinpsi;
    Rot(1,0) = sinphi*costheta*cospsi + cosphi*sinpsi;
    Rot(2,0) = -sintheta*cospsi;

    Rot(0,1) = -cosphi*costheta*sinpsi - sinphi*cospsi;
    Rot(1,1) = -sinphi*costheta*sinpsi + cosphi*cospsi;
    Rot(2,1) = sintheta*sinpsi;

    Rot(0,2) = cosphi*sintheta;
    Rot(1,2) = sinphi*sintheta;
    Rot(2,2) = costheta;
    return (Rot);
  }

// Euler rotation matrix
// ZYZ convention - passive rotation
// to rotate a colvec V, do Rt * V
// [[Rcpp::export]]
arma::mat cpp_euler_passive(const double phi, const double theta, const double psi)
  {
    arma::mat Rot(3,3);
    const double cosphi = cos(phi), cospsi = cos(psi), costheta = cos(theta);
    const double sinphi = sin(phi), sinpsi = sin(psi), sintheta = sin(theta);
    Rot(0,0) = cosphi*costheta*cospsi - sinphi*sinpsi;
    Rot(0,1) = sinphi*costheta*cospsi + cosphi*sinpsi;
    Rot(0,2) = -sintheta*cospsi;

    Rot(1,0) = -cosphi*costheta*sinpsi - sinphi*cospsi;
    Rot(1,1) = -sinphi*costheta*sinpsi + cosphi*cospsi;
    Rot(1,2) = sintheta*sinpsi;

    Rot(2,0) = cosphi*sintheta;
    Rot(2,1) = sinphi*sintheta;
    Rot(2,2) = costheta;
    return (Rot);
  }

// Rotation matrix about a cartesian axis
// [[Rcpp::export]]
arma::mat cpp_axis_rotation(const double angle, const int axis)
  {
    arma::mat Rot(3,3);
    const double ca = cos(angle), sa = sin(angle);

    if(axis == 0) {// rotate about x axis
      Rot(0,0) = 1;   Rot(0,1) = 0 ;  Rot(0,2) =   0;
      Rot(1,0) = 0;   Rot(1,1) = ca;  Rot(1,2) = -sa;
      Rot(2,0) = 0;   Rot(2,1) = sa;  Rot(2,2) =  ca;
    }
    if(axis == 1){// rotate about y axis
      Rot(0,0) =  ca; Rot(0,1) = 0;   Rot(0,2) = sa;
      Rot(1,0) =  0 ; Rot(1,1) = 1;   Rot(1,2) =  0;
      Rot(2,0) = -sa; Rot(2,1) = 0;   Rot(2,2) = ca;
    }
    if(axis == 2){// rotate about z axis
      Rot(0,0) = ca;  Rot(0,1) = -sa; Rot(0,2) = 0;
      Rot(1,0) = sa;  Rot(1,1) =  ca; Rot(1,2) = 0;
      Rot(2,0) =  0;  Rot(2,1) =   0; Rot(2,2) = 1;
    }

    return (Rot);
  }
