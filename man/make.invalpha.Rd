\name{make.invalpha}
\alias{make.invalpha}
\title{make.invalpha}
\usage{
  make.invalpha(cluster, material,
    polarizability.fun = polarizability.ellipsoid, ...)
}
\arguments{
  \item{cluster}{cluster}

  \item{material}{material}

  \item{polarizability.fun}{polarizability function}

  \item{...}{additional arguments passed to
  polarizability.fun}
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

