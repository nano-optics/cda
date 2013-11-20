\name{cluster_dimer2}
\alias{cluster_dimer2}
\title{cluster_dimer2}
\usage{
cluster_dimer2(d = a, dihedral = 0, alpha1 = 0, alpha2 = 0, a = 0.035,
  b = 0.012)
}
\arguments{
  \item{d}{center-to-center distance}

  \item{dihedral}{dihedral angle}

  \item{alpha1}{angle first rod}

  \item{alpha2}{angle second rod}

  \item{a}{semi axis}

  \item{b}{semi axis}
}
\value{
  list with r, sizes, angles
}
\description{
  cluster_dimer2
}
\details{
  cluster with two nanorods first rod along x at (0, 0,
  -d/2) second rod at (0, 0, d/2)
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

