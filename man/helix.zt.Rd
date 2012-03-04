\name{helix.zt}
\alias{helix.zt}
\title{helix2}
\usage{
  helix.zt(R0 = 0.5, pitch = 0.6,
    z = runif(100, -0.5, 0.5),
    theta = runif(100, 0, 2 * pi), right = TRUE)
}
\arguments{
  \item{R0}{radius in um}

  \item{pitch}{pitch in um}

  \item{z}{positions of particles along cylinder axis}

  \item{theta}{angular position along cylinder}

  \item{right}{logical, handedness}
}
\value{
  matrix N x 3
}
\description{
  helix.zt
}
\details{
  create a matrix of xyz positions for ellipsoids arranged
  on a cylinder with an helical twist
}
\author{
  baptiste Auguie
}

