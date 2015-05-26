library(cda)

wavelength <- 500
medium <- 1.33
kn <- 2*pi/wavelength*medium

E0L=1/sqrt(2)*c(0,1,1i)
E0R=1/sqrt(2)*c(0,1i,1)
k0=c(1,0,0)

cl <- cluster_chain(5000)
Beta <- inverse_polarizability(cl, material=epsAu(wavelength), medium=medium)

Angles <- cbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2

A <- cda$interaction_matrix(cl$r, kn, c(Beta), cl$angles, TRUE)
Ei <- cda$incident_field(E0L, k=kn*k0, r=cl$r, Angles)

system.time(P2 <- cda$cg_solve(A, Ei,  0*Ei, 10, 1e-5))
# P <- solve(A, Ei)
# range(abs(P-P2)/max(Mod(P)))

library(microbenchmark)

microbenchmark(r = solve(A, Ei),
               cpp = cda$cg_solve(A, Ei,  0*Ei, 10, 1e-5), times = 5)

# for 10 dipoles
# Unit: microseconds
# expr     min      lq      mean  median      uq     max neval
# r  30.055  33.960  36.47319  34.539  36.641 138.704   100
# cpp 128.518 130.757 141.13337 138.051 144.046 253.690   100

# for 50 dipoles
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# r 2.167841 2.283330 2.523692 2.433572 2.731529 3.467672   100
# cpp 1.575142 1.763508 2.027794 1.995655 2.200366 4.733692   100

# for 100 dipoles
# Unit: milliseconds
# expr       min        lq      mean    median        uq       max neval
# r 16.603467 17.073234 18.394929 18.218043 19.459713 23.096228   100
# cpp  6.381693  6.868462  7.705662  7.657668  8.314182  9.937314   100

# for 1000 dipoles
# Unit: milliseconds
# expr        min         lq       mean     median         uq       max neval
# r 16993.1205 17979.3306 18921.1487 18842.4223 20125.4297 20504.922    10
# cpp   815.6591   865.8921   904.8368   880.1814   968.7177  1010.624    10