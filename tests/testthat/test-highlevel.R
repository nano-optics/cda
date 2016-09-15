library(cda)
library(testthat)

context("Checking low-level functions")

material <- data.frame(wavelength=500, epsilon=-10+0.2i)
medium <- 1.0
kn <- 2*pi/material$wavelength * medium
kvec = c(0,0,1)
Evec <- c(1,0,0)
cl <- list(positions = cbind(c(0, 0, 0),
                             c(100, 0, 0)),
           angles = cbind(c(pi/3, 0, 0),
                          c(0, 0, 0)),
           sizes = cbind(c(30, 10, 10),
                         c(30, 10, 10)))

Alpha <- alpha_ellipsoid(cl$sizes, material, medium)
Incidence <- matrix(0, ncol=3,nrow=4)

AlphaBlocks <- cda:::cpp_alpha_blocks(Alpha, cl$angles)

A <- cda:::cpp_interaction_matrix(cl$positions, as.double(kn), AlphaBlocks)

Ein <- cda:::cpp_incident_field(Evec, kvec*kn, cl$positions, Incidence)

E <- solve(A, Ein)
P <- cda:::cpp_polarization(E, AlphaBlocks)

cext <- cda:::cpp_extinction(kn, P, Ein)
cabs <- cda:::cpp_absorption(kn, P, E)

Scattering <- quadrature_sphere(Nq=50, "gl")

res1 <- cda:::cpp_dispersion_spectrum(kn, medium, cl$positions, Alpha, 
                                      cl$angles, matrix(0.0), 0L, 
                                      Scattering$nodes, Scattering$weights, 
                                      0L, 2L, 
                                      50, 1e-4, FALSE)


test_that("High-level dispersion spectrum coincides with low-level",
          expect_equal(res1[1:3], c(cext[1], cabs[1], cext[1] - cabs[1])/2, tolerance = 1e-3)
)

