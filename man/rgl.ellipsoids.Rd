\name{rgl.ellipsoids}
\alias{rgl.ellipsoids}
\title{rgl.ellipsoids}
\usage{
  rgl.ellipsoids(positions, sizes, angles, ...)
}
\arguments{
  \item{positions}{matrix of positions}

  \item{sizes}{matrix of axis lengths}

  \item{angles}{matrix of Euler angles}

  \item{...}{additional params}
}
\value{
  rgl mesh
}
\description{
  Create a list of rgl ellipsoids oriented in space
}
\details{
  each ellipsoid is specified by its position, dimensions,
  and Euler angles
}
\examples{
cl <- helix(0.5, 1, 36, delta=pi/6, n.smooth=1e3)
sizes <- clust.equalsizes(0.04,0.02,0.02,NROW(cl$positions))
\dontrun{ require(rgl) ; rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold") }
}
\author{
  baptiste Auguie
}

