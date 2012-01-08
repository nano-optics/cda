#ifndef _cda_UTILS_H
#define _cda_UTILS_H

#include <RcppArmadillo.h>

arma::mat euler(const double phi, const double theta, const double psi);
void progress_bar(double x, double N);

#endif
