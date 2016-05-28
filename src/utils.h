#ifndef _UTILS_H
#define _UTILS_H

arma::mat cpp_euler_active(const double phi, const double theta, const double psi);
arma::mat cpp_euler_passive(const double phi, const double theta, const double psi);
arma::mat cpp_axis_rotation(const double angle, const int axis);
void progress_bar(double x, double N);

#endif
