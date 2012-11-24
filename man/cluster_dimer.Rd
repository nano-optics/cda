\name{cluster_dimer}
\alias{cluster_dimer}
\title{cluster_dimer}
\usage{
  cluster_dimer(d = a, dihedral = 0, alpha1 = 0,
    alpha2 = 0, a = 0.035, b = 0.012)
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
  cluster_dimer
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
  \code{\link{cluster_helix}}, \code{\link{equal_angles}},
  \code{\link{equal_sizes}}, \code{\link{helix}}
}

