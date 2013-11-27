\name{inverse_polarizability}
\alias{inverse_polarizability}
\title{inverse_polarizability}
\usage{
inverse_polarizability(cluster, material,
  polarizability_fun = polarizability_ellipsoid, ...)
}
\arguments{
  \item{cluster}{cluster}

  \item{material}{material}

  \item{polarizability_fun}{polarizability function}

  \item{...}{additional arguments passed to
  polarizability_fun}
}
\value{
matrix with each row being the 3 principal values of each
polarizability tensor
}
\description{
inverse polarizability tensors
}
\details{
calculates and formats the principal polarizability of
several particles
}
\author{
Baptiste Auguie
}
\seealso{
Other user_level polarizability:
\code{\link{alpha_kuwata}}, \code{\link{La}},
\code{\link{polarizability_ellipsoid}}
}

