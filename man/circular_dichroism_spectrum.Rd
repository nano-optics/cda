\name{circular_dichroism_spectrum}
\alias{circular_dichroism_spectrum}
\title{circular_dichroism_spectrum}
\usage{circular_dichroism_spectrum(cluster, material, n=1.33, N=100, averaging=c("QMC","GL","grid"), full=TRUE, progress=FALSE,
    result.matrix=FALSE)
}
\description{Simulate a CD spectrum}
\details{CD spectrum}
\value{data.frame or matrix of CD spectrum}
\author{baptiste Auguie}
\arguments{\item{cluster}{cluster (list)}
\item{material}{material}
\item{n}{refractive index medium}
\item{N}{number of integration points}
\item{averaging}{averaging method, using either Quasi Monte Carlo, Gauss Legendre, or regular grid}
\item{full}{logical use full (retarded) dipolar field}
\item{progress}{print progress lines}
\item{result.matrix}{logical return the results as a matrix}
}

