\name{cda-module}
\title{
	Rcpp module: cda
}
\description{
	Exposes basic C++ functions used in the coupled-dipole approximation.
}
\details{
	The module contains the following items: 
	
	functions: \describe{
        \item{absorption}{absorption cross-section}
        \item{euler}{ 3D rotation matrix parametrized by Euler angles }
        \item{extinction}{ extinction cross-section }
        \item{interaction_matrix}{ build the coupled-dipole interaction matrix }
		} 
	
	
}
\examples{
show( cda )
}
\keyword{datasets}
