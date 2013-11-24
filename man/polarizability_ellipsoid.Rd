\name{polarizability_ellipsoid}
\alias{polarizability_ellipsoid}
\title{polarizability_ellipsoid}
\usage{
polarizability_ellipsoid(wavelength, epsilon, a = 50, b = 30, c = b,
  medium = 1.33, kuwata = TRUE)
}
\arguments{
  \item{wavelength}{wavelength in nm}

  \item{epsilon}{complex permittivity}

  \item{a}{semi-axis in nm}

  \item{b}{semi-axis in nm}

  \item{c}{semi-axis in nm}

  \item{medium}{surrounding medium}

  \item{kuwata}{logical, use Kuwata or Clausius Mossotti
  prescription, see Details}
}
\value{
matrix of polarizability
}
\description{
principal polarizability components for an ellipsoidal
particle
}
\details{
uses the Kuwata prescription (see references)

The Kuwata version includes semi-empirical terms of
radiative correction and dynamic depolarisation to better
match the fully retarded dipolar response in a reasonable
range of (subwavelength) sizes and aspect ratios.
}
\author{
baptiste Auguie
}
\references{
Kuwata et al. Resonant light scattering from metal
nanoparticles: Practical analysis beyond Rayleigh
approximation Appl. Phys. Lett. 83, 22 (2003)
}
\seealso{
Other user_level polarizability:
\code{\link{alpha_kuwata}},
\code{\link{inverse_polarizability}}, \code{\link{La}}
}

