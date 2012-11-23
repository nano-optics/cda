\name{alpha_kuwata}
\alias{alpha_kuwata}
\alias{Kuwata.A}
\alias{Kuwata.B}
\title{alpha_kuwata}
\usage{
  alpha_kuwata(wavelength, epsilon, V, axis, L,
    medium = 1.33)
}
\arguments{
  \item{wavelength}{wavelength}

  \item{epsilon}{permittivity}

  \item{V}{volume}

  \item{axis}{semi-axis along incident field}

  \item{L}{shape factor}

  \item{medium}{refractive index}
}
\value{
  polarizability
}
\description{
  polarizability
}
\details{
  prescription from Kuwata
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level polarizability:
  \code{\link{inverse_polarizability}}, \code{\link{La}},
  \code{\link{polarizability_ellipsoid}}
}

