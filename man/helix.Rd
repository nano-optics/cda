\name{helix}
\alias{helix}
\title{helix}
\usage{
  helix(R0 = 0.5, pitch = 0.6, N = 5, delta = pi/8,
    delta0 = pi/2, n.smooth = 100 * N, right = TRUE)
}
\arguments{
  \item{R0}{radius in microns}

  \item{pitch}{pitch in microns}

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
cl <- helix(0.5, 1, 36, delta=pi/6, n.smooth=1e3) ; str(cl)
\dontrun{require(rgl)
open3d()
spheres3d(cl$smooth, radius=0.01,col=2)
spheres3d(cl$positions, radius=0.1, col="gold")
## ellipsoids are oriented following the helix
sizes <- clust.equalsizes(0.04,0.02,0.02,NROW(cl$positions))
rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold") }
}
\author{
  baptiste Auguie
}

