\name{dispersion_spectrum}
\alias{dispersion_spectrum}
\title{dispersion_spectrum}
\usage{
  dispersion_spectrum(cluster, angles, axis = "z",
    material, medium = 1.33,
    polarisation = c("linear", "circular"), invert = FALSE,
    progress = FALSE)
}
\arguments{
  \item{cluster}{list describing a cluster}

  \item{angles}{of incident field in radians}

  \item{axis}{of incident field rotation character vector
  from ('x', 'y', 'z')}

  \item{material}{list}

  \item{medium}{medium refractive index}

  \item{polarisation}{linear or circular polarisation}

  \item{invert}{logical}

  \item{progress}{logical, display progress bar}
}
\value{
  data.frame
}
\description{
  dispersion spectrum
}
\details{
  dispersion spectrum
}
\author{
  baptiste Auguie
}

