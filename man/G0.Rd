\name{G0}
\alias{G0}
\alias{gfun}
\docType{data}
\title{Lattice sum}
\description{Converged lattice sum G0 for an infinite 2D array of dipoles at normal incidence.}
\usage{data(G0)}
\format{
  A data frame with 1000 observations on the following 3 variables.
  \describe{
    \item{\code{wavelength}}{a numeric vector}
    \item{\code{Qx}}{a numeric vector}
    \item{\code{Gxx}}{a complex vector}
  }
}
\details{The calculation was made using code from Prof. J. G. de Abajo (CSIC, Spain)}
\examples{
data(G0)
\dontrun{demo(lattice_sum)}
}
\keyword{datasets}
