\docType{data}
\name{cda}
\alias{cda}
\title{Rcpp module: cda}
\description{
Exposes basic C++ functions used in the coupled-dipole
approximation.
}
\details{
\itemize{ \item{absorption}{ Absorption cross-section}
\item{extinction}{ Extinction cross-section }
\item{axis_rotation}{ 3D rotation matrix parametrized by
axis + angle} \item{euler}{ 3D rotation matrix parametrized
by Euler angles } \item{interaction_matrix}{ Build the
coupled-dipole interaction matrix } \item{block_diagonal}{
Diagonal blocks of the coupled-dipole interaction matrix }
\item{incident_field}{ Construct the incident fields for
specific Euler angles} \item{multiple_incident_field}{
Construct the incident fields for specific axes+angles} }
}
\examples{
show( cda )
}

