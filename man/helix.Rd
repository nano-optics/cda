\name{helix}
\alias{helix}
\title{helix}
\usage{
  helix(R0 = 500, pitch = 600, N = 5, delta = pi/8,
    delta0 = pi/2, n.smooth = 100 * N, right = TRUE)
}
\arguments{
  \item{R0}{radius in nm}

  \item{pitch}{pitch in nm}

  \item{N}{number of particles}

  \item{delta}{twist angle between two particles}

  \item{delta0}{angle shift at origin}

  \item{n.smooth}{number of interpolation points (for
  plotting purposes)}

  \item{right}{logical, handedness}
}
\value{
  list with r, sizes, invalpha, angles, R0 and smooth
  interp. points
}
\description{
  helix curve
}
\details{
  add particles on an helix
}
\examples{
cl <- helix(500, 1000, 36, delta=pi/6, n.smooth=1e3) ; str(cl)
\dontrun{
require(rgl)
open3d()
spheres3d(cl$smooth, radius=1,col=2)
## ellipsoids are oriented following the helix
sizes <- equal_sizes(40, 20, 20,NROW(cl$positions))
rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold")
}
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level cluster: \code{\link{cluster_chain}},
  \code{\link{cluster_dimer}}, \code{\link{cluster_helix}},
  \code{\link{equal_angles}}, \code{\link{equal_sizes}}
}

