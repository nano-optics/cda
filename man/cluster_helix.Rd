\name{cluster_helix}
\alias{cluster_helix}
\title{cluster_helix}
\usage{
  cluster_helix(N = 5, R0 = 12, pitch = 15, delta = pi/2,
    delta0 = 0, right = TRUE, a = 5, b = a, c = a,
    angles = c("helix", "random", "fixed"), seed = 123,
    ...)
}
\arguments{
  \item{N}{number of particles}

  \item{R0}{radius of helix}

  \item{pitch}{pitch of helix}

  \item{delta}{angle between particles}

  \item{delta0}{initial angle}

  \item{right}{logical, helicity}

  \item{a}{ellipsoid semi-axis}

  \item{b}{ellipsoid semi-axis}

  \item{c}{ellipsoid semi-axis}

  \item{angles}{type of angular orientation}

  \item{seed}{random seed for reproducibility}

  \item{...}{extra arguments (ignored)}
}
\value{
  list
}
\description{
  cluster_helix
}
\details{
  helical cluster of ellipsoids
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level cluster: \code{\link{cluster_chain}},
  \code{\link{cluster_dimer}}, \code{\link{equal_angles}},
  \code{\link{equal_sizes}}, \code{\link{helix}}
}

