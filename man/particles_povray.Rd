\name{particles_povray}
\alias{particles_povray}
\title{particles_povray}
\usage{
  particles_povray(positions, angles, sizes,
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
  Other user_level povray: \code{\link{curve_povray}}
}

