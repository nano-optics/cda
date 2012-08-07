\name{makeSpheresCluster}
\alias{makeSpheresCluster}
\title{makeSpheresCluster}
\usage{
  makeSpheresCluster(radius = 0.005, N, R0 = 0.5,
    pitch = 0.6, delta = pi/8, delta0 = 0, right = TRUE)
}
\arguments{
  \item{radius}{radius}

  \item{N}{number of spheres}

  \item{R0}{radius helix}

  \item{pitch}{pitch helix}

  \item{delta}{phase helix}

  \item{delta0}{initial phase helix}

  \item{right}{logical, handedness}
}
\value{
  list of positions, sizes and angles
}
\description{
  make a cluster of spheres
}
\details{
  make a cluster of spheres
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level cluster:
  \code{\link{clust.equalangles}},
  \code{\link{clust.equalsizes}}, \code{\link{helix}},
  \code{\link{helix.zt}}, \code{\link{makeDimerCluster}},
  \code{\link{makeDimerDihedral}},
  \code{\link{makeHelixCluster}},
  \code{\link{makeHelixCluster.zt}},
  \code{\link{makeRodChain}}
}

