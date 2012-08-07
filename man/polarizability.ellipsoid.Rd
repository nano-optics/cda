\name{polarizability.ellipsoid}
\alias{polarizability.ellipsoid}
\title{polarizability.ellipsoid}
\usage{
  polarizability.ellipsoid(lambda, epsilon, a = 0.05,
    b = 0.03, c = b, n = 1.33, kuwata = TRUE)
}
\arguments{
  \item{lambda}{wavelength in microns}

  \item{epsilon}{complex permittivity}

  \item{a}{semi-axis in um}

  \item{b}{semi-axis in um}

  \item{c}{semi-axis in um}

  \item{n}{surrounding medium}

  \item{kuwata}{logical, use Kuwata or Clausius Mossotti}
}
\value{
  matrix of polarizability
}
\description{
  principal polarizability components for an ellipsoidal
  particle
}
\details{
  uses the Kuwata prescription
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level polarizability: \code{\link{La}},
  \code{\link{alpha_kuwata}}, \code{\link{make.invalpha}}
}

