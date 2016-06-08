#ifndef _INCIDENT_H
#define _INCIDENT_H

arma::cx_mat cpp_incident_field(const arma::cx_colvec& E0,
                                const arma::colvec& k,
                                const arma::mat& R,
                                const arma::mat& Angles);

arma::cx_mat cpp_incident_field_axis(const arma::cx_colvec& Evec,
                                     const arma::colvec& kvec,
                                     const arma::mat& R,
                                     const arma::colvec& Incidence,
                                     const arma::ivec& Axes);

#endif
