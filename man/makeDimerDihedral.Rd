\name{makeDimerDihedral}
\alias{makeDimerDihedral}
\title{makeDimerDihedral}
\usage{
  makeDimerDihedral(d = a, dihedral = 0, alpha1 = 0,
    alpha2 = 0, a = 0.035, b = 0.012, right = TRUE)
}
\arguments{
  \item{d}{center-to-center distance}

  \item{dihedral}{dihedral angle}

  \item{alpha1}{angle first rod}

  \item{alpha2}{angle second rod}

  \item{a}{semi axis}

  \item{b}{semi axis}

  \item{right}{handedness}
}
\value{
  list with r, sizes, angles
}
\description{
  makeDimerDihedral
}
\details{
  makeDimerDihedral first rod along x at (0, 0, -d/2)
  second rod at (0, 0, d/2)
}
\author{
  baptiste Auguie
}

