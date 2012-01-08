\name{rgl.ellipsoid}
\alias{rgl.ellipsoid}
\title{rgl.ellipsoid}
\usage{rgl.ellipsoid(x=0, y=0, z=0, a=1, b=1, c=1, phi=0, theta=0, psi=0, subdivide=3,
    smooth=TRUE, ...)
}
\description{creates an rgl ellipsoid}
\details{deforms, rotate, and translate a sphere}
\value{an rgl mesh}
\author{baptiste Auguie}
\arguments{\item{x}{x}
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
\examples{\dontrun{ require(rgl) ;  ee <- rgl.ellipsoid()
shapelist3d(ee) }}

