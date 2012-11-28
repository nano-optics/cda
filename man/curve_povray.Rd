\name{curve_povray}
\alias{curve_povray}
\title{curve_povray}
\usage{
  curve_povray(positions, size = 0.005, radius,
    out = "positions.pov")
}
\arguments{
  \item{positions}{matrix of positions}

  \item{size}{radius of particles}

  \item{radius}{radius of inner cylinder}

  \item{out}{output filename}
}
\value{
  side-effect only (note append=TRUE)
}
\description{
  create long string of small particles for povray input
}
\details{
  writes a list of particles, as well as cylinder
}
\author{
  baptiste Auguie
}
\seealso{
  Other user_level povray: \code{\link{particles_povray}}
}

