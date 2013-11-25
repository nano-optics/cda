\name{circular_dichroism_spectrum}
\alias{circular_dichroism_spectrum}
\title{circular_dichroism_spectrum}
\usage{
circular_dichroism_spectrum(cluster, material, medium = 1.33, Nquad = 100,
  averaging = c("GL", "QMC", "grid", "cheap"), iterative = FALSE,
  precision = 0.001, Nmax = 10000, dN = Nquad, full = TRUE,
  progress = FALSE, verbose = TRUE, result.matrix = FALSE)
}
\arguments{
  \item{cluster}{cluster (list)}

  \item{material}{material}

  \item{medium}{refractive index medium}

  \item{Nquad}{number of integration points}

  \item{averaging}{averaging method, using either Gauss
  Legendre quadrature (default), Quasi Monte Carlo, regular
  grid, or "cheap" (3 axes)}

  \item{iterative}{logical, increase N until convergence
  (QMC only)}

  \item{precision}{relative diff between two runs (QMC
  only)}

  \item{Nmax}{maximum N if convergence not attained (QMC
  only)}

  \item{dN}{iterative increase in N (QMC only)}

  \item{full}{logical use full (retarded) dipolar field}

  \item{progress}{print progress lines}

  \item{verbose}{display messages}

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
\references{
Y. Okada, Efficient numerical orientation averaging of
light scattering properties with a quasi-Monte-Carlo
method, Journal of Quantitative Spectroscopy and Radiative
Transfer, Volume 109, Issue 9, June 2008, Pages 1719-1742.
}

