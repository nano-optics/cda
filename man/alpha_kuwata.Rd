\name{alpha_kuwata}
\alias{Kuwata.A}
\alias{Kuwata.B}
\alias{alpha_kuwata}
\title{alpha_kuwata}
\usage{
  alpha_kuwata(lambda, epsilon, V, axis, L,
    epsilon.medium = 1.33^2)
}
\arguments{
  \item{lambda}{wavelength}

  \item{epsilon}{permittivity}

  \item{V}{volume}

  \item{axis}{semi-axis along incident field}

  \item{L}{shape factor}

  \item{epsilon.medium}{medium permittivity}
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
  Other user_level polarizability: \code{\link{La}},
  \code{\link{make.invalpha}},
  \code{\link{polarizability.ellipsoid}}
}

