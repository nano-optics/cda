\name{rgl.ellipsoid}
\alias{rgl.ellipsoid}
\title{rgl.ellipsoid}
\usage{
rgl.ellipsoid(x = 0, y = 0, z = 0, a = 1, b = 1, c = 1, phi = 0,
  theta = 0, psi = 0, subdivide = 3, smooth = TRUE, ...)
}
\arguments{
  \item{x}{x}

  \item{y}{y}

  \item{z}{z}

  \item{a}{axis}

  \item{b}{axis}

  \item{c}{axis}

  \item{phi}{phi}

  \item{theta}{theta}

  \item{psi}{psi}

  \item{subdivide}{subdivision}

  \item{smooth}{smoothing}

  \item{...}{additional params}
}
\value{
an rgl mesh
}
\description{
creates an rgl ellipsoid
}
\details{
deforms, rotate, and translate a sphere
}
\examples{
\dontrun{ require(rgl) ;  ee <- rgl.ellipsoid()
shapelist3d(ee) }
}
\author{
baptiste Auguie
}
\seealso{
Other user_level rgl: \code{\link{rgl.ellipsoids}}
}

