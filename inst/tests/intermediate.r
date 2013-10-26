library(testthat)
library(cda)
context("Checking the intermediate quantities in CD calculations")

wavelength <- c(500, 600)
material <- epsAu(wavelength)
medium <- 1.3  

k0 <- 2*pi/wavelength
kn <- k0*medium

cluster <- list(r = rbind(c(0, 0, 0),
                     c(0, 0, 200)),
           angles = rbind(c(0, 0, 0),
                          c(pi/4, 0, 0)),
           sizes = rbind(c(40, 20, 20),
                         c(40, 20, 20)))

invalpha <- inverse_polarizability(cluster, material, 
                                  polarizability_fun=polarizability_ellipsoid, 
                                  medium=medium, kuwata=TRUE)

.invalpha <- structure(c(-1.82237397053259e-05-3.56779387878792e-05i, -2.39972560681825e-06-5.46975933994133e-06i, 
                         2.89922059965679e-05-3.56779387878792e-05i, 4.41715911001179e-05-5.4697593399413e-06i, 
                         2.89922059965679e-05-3.56779387878792e-05i, 4.41715911001179e-05-5.4697593399413e-06i, 
                         -1.82237397053259e-05-3.56779387878792e-05i, -2.39972560681825e-06-5.46975933994133e-06i, 
                         2.89922059965679e-05-3.56779387878792e-05i, 4.41715911001179e-05-5.4697593399413e-06i, 
                         2.89922059965679e-05-3.56779387878792e-05i, 4.41715911001179e-05-5.4697593399413e-06i
), .Dim = c(2L, 6L), .Dimnames = list(NULL, c("alpha.kuwata.a", 
                                              "alpha.kuwata.b", "alpha.kuwata.c", "alpha.kuwata.a", "alpha.kuwata.b", 
                                              "alpha.kuwata.c")))

test_that("the inverse polarisability is the same as earlier versions", {
  expect_equal(invalpha,.invalpha)
})

## cheap averaging

nodes <- rbind(c(1/2, 0), # +x is phi=0, psi=0
               c(1/2, 1/4), # +y is phi=pi/2, psi=0
               c(1, 1/4)) # +z is phi=pi/2, psi=pi/2

res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, 
                                       as.matrix(nodes),
                                       as.integer(TRUE), as.integer(TRUE))

.res <- structure(c(7629.71366621336, 18093.6538420437, 6965.12045086035, 
                   12499.019321288, -0.216853598901253, -533.349596245196, 2.67479400633874, 
                   -361.940147617599), .Dim = c(2L, 4L))


test_that("the result of cheap orientation averaging is the same as earlier versions", {
  expect_equal(res,.res)
})


Beta <- matrix(invalpha[1,], nrow=3)

A <- cda$interaction_matrix(cluster$r, kn[1], Beta, cluster$angles, 1)
.A <- structure(c(-1.82237397053259e-05-3.56779387878792e-05i, 0+0i, 
                  0+0i, 1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 
                  0+0i, 2.89922059965679e-05-3.56779387878792e-05i, 0+0i, 0+0i, 
                  1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 0+0i, 
                  2.89922059965679e-05-3.56779387878792e-05i, 0+0i, 0+0i, 3.5040262644085e-07-7.79039958472603e-07i, 
                  1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 5.384233145621e-06-3.56779387878792e-05i, 
                  -0.0000236079728509469+0i, 0+0i, 0+0i, 1.14864729173871e-06+5.5676095065754e-07i, 
                  0+0i, -0.0000236079728509469+0i, 5.384233145621e-06-3.56779387878792e-05i, 
                  0+0i, 0+0i, 0+0i, 3.5040262644085e-07-7.79039958472603e-07i, 
                  0+0i, 0+0i, 2.89922059965679e-05-3.56779387878792e-05i), .Dim = c(6L, 
                                                                                    6L))

all.equal(A, .A)
max(Re(A))
