\name{circular_dichroism_spectrum}
\alias{circular_dichroism_spectrum}
\title{circular_dichroism_spectrum}
\usage{
  circular_dichroism_spectrum(cluster, material,
    medium = 1.33, N = 100,
    averaging = c("GL", "QMC", "grid"), full = TRUE,
    progress = FALSE, result.matrix = FALSE)
}
\arguments{
  \item{cluster}{cluster (list)}

  \item{material}{material}

  \item{medium}{refractive index medium}

  \item{N}{number of integration points}

  \item{averaging}{averaging method, using either Gauss
  Legendre quadrature (default), Quasi Monte Carlo, or
  regular grid}

  \item{full}{logical use full (retarded) dipolar field}

  \item{progress}{print progress lines}

  \item{result.matrix}{logical return the results as a
  matrix}
}
\description{
  Simulate a CD spectrum
}
\details{
  CD spectrum
}
\author{
  baptiste Auguie
}

