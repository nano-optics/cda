\name{array_factor}
\alias{array_factor}
\title{array factor}
\usage{
  array_factor(wavelength, N, pitch)
}
\arguments{
  \item{wavelength}{wavelength in nm}

  \item{N}{half number of dipoles along one side}

  \item{pitch}{pitch in nm}
}
\value{
  S
}
\description{
  C++ calculation of the array factor
}
\details{
  C++ calculation of the array factor
}
\examples{
S <- array_factor(seq(400, 600),  10,  500)
str(S)
}
\author{
  baptiste Auguie
}

