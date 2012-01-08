\name{makeDimerCluster}
\alias{makeDimerCluster}
\title{makeDimerCluster}
\usage{makeDimerCluster(d=a, phi=pi/4, psi=0, a=0.035, b=0.012, c=b, right=TRUE)
}
\description{makeDimerCluster}
\details{makeDimerCluster
first rod along x at (0, 0, -d/2)
second rod at (0, 0, d/2),  rotated by theta (z) in [0, pi], psi (z') in [-pi/2, pi/2]}
\value{list with r,  sizes,  angles}
\author{baptiste Auguie}
\arguments{\item{d}{center-to-center distance}
\item{phi}{longitude}
\item{psi}{latitude}
\item{a}{semi axis}
\item{b}{semi axis}
\item{c}{semi axis}
\item{right}{logical,  handedness}
}

