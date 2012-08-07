\name{particles.povray}
\alias{particles.povray}
\title{particles.povray}
\usage{
  particles.povray(positions, angles, sizes,
    out = "positions.pov")
}
\arguments{
  \item{positions}{matrix of positions}

  \item{angles}{matrix of Euler angles in radians}

  \item{sizes}{matrix of particle sizes}

  \item{out}{output filename}
}
\value{
  side-effect only
}
\description{
  create particles for povray input
}
\details{
  writes a list of particles
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level povray: \code{\link{curve.povray}}
}

