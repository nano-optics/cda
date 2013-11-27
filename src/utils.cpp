#include "utils.h"
#include <RcppArmadillo.h>
#include <iostream>

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
arma::mat euler(const double phi, const double theta, const double psi)
  {
    arma::mat Rot(3,3);
    const double cosphi = cos(phi), cospsi = cos(psi), costheta = cos(theta);
    const double sinphi = sin(phi), sinpsi = sin(psi), sintheta = sin(theta);
    Rot(0,0) = cospsi*cosphi - costheta*sinphi*sinpsi;
    Rot(0,1) = cospsi*sinphi + costheta*cosphi*sinpsi; 
    Rot(0,2) = sinpsi*sintheta;
       
    Rot(1,0) = -sinpsi*cosphi - costheta*sinphi*cospsi; 
    Rot(1,1) = -sinpsi*sinphi + costheta*cosphi*cospsi; 
    Rot(1,2) = cospsi*sintheta;
       
    Rot(2,0) = sinphi*sintheta;
    Rot(2,1) = -cosphi*sintheta; 
    Rot(2,2) = costheta;
    return (Rot);
  }

// Rotation matrix about a cartesian axis
 arma::mat axis_rotation(const double angle, const int axis)
  {
    arma::mat Rot(3,3);
    const double cosangle = cos(angle), sinangle = sin(angle);

    if(axis == 0) {// rotate about x axis
      Rot << 1 << 0 << 0 << arma::endr
	<< 0 << cosangle << -sinangle << arma::endr
	<< 0 << sinangle << cosangle << arma::endr;

    }
    if(axis == 1){// rotate about y axis
    Rot << cosangle << 0 << sinangle << arma::endr
	<< 0 << 1 << 0 << arma::endr
	<< -sinangle << 0 << cosangle << arma::endr;

    }
    if(axis == 2){// rotate about z axis
    Rot << cosangle << -sinangle << 0 << arma::endr
	<< sinangle << cosangle << 0 << arma::endr
	<< 0 << 0 << 1 << arma::endr;

    }

    return (Rot);
  }


