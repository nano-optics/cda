\name{cluster_dimer_end}
\alias{cluster_dimer_end}
\title{cluster_dimer_end}
\usage{
cluster_dimer_end(d = a, dihedral = 0, a = 0.035, b = 0.012,
  rescale = TRUE)
}
\arguments{
  \item{d}{end-to-end distance}

  \item{dihedral}{dihedral angle}

  \item{a}{semi axis}

  \item{b}{semi axis}

  \item{rescale}{logical, rescale the z coordinates so that
  d is the center-to-center distance}
}
\value{
list with r, sizes, angles
}
\description{
cluster_dimer_end
}
\details{
cluster with two nanorods first rod along x at (0, 0, -d/2)
second rod at (0, 0, d/2)
}
\author{
baptiste Auguie
}
\seealso{
Other user_level cluster: \code{\link{cluster_chain}},
\code{\link{cluster_dimer}}, \code{\link{cluster_helix}},
\code{\link{equal_angles}}, \code{\link{equal_sizes}},
\code{\link{helix}}
}

