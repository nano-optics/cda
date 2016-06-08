#ifndef _CROSSSECTIONS_H
#define _CROSSSECTIONS_H

arma::colvec cpp_extinction(const double kn, const arma::cx_mat& P,
                            const arma::cx_mat& E0);

arma::colvec cpp_absorption(const double kn, const arma::cx_mat& P,
                            const arma::cx_mat& E);

arma::colvec cpp_scattering(const arma::mat& positions,
                            const arma::mat& angles,
                            const arma::colvec& weights,
                            const double kn,
                            const arma::cx_mat& P);

#endif
