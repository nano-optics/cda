#setwd("~/Documents/plasmonics/cd/vignettes/userguide")

library(cd)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)


# dielectric function
# wvl <- seq(500, 750, length=200)
# gold <- epsAu(wvl)
# cl <- cluster_dimer(100)

# Rprof("profiling.out")
material <- data.frame(wavelength=500, epsilon=-10+0.2i)
medium <- 1.0
kn <- 2*pi/material$wavelength * medium
kvec = c(0,0,1)
Evec <- c(1,0,0)
N <- 1000
cl <- cluster_chain(N)

Alpha <- alpha_ellipsoid(cl$sizes, material, medium)
AlphaBlocks <- cd:::cpp_alpha_blocks(Alpha, cl$angles)

R <- cl$positions
full <- TRUE

A <- cd:::cpp_interaction_matrix(R, as.double(kn), AlphaBlocks, full)

Incidence <- matrix(0, ncol=3,nrow=4)

Ein <- cd:::cpp_incident_field(Evec, kvec*kn, R, Incidence)

# switch to accelerate
# https://groups.google.com/d/msg/r-sig-mac/k4rDRRdtNwE/oIvxyd0Yd88J
# http://blog.quadrivio.com/2015/06/improved-r-performance-with-openblas.html
system.time(E <- solve(A, Ein))
# user  system elapsed 
# 17.042   0.098  17.187 
# vecLib
# user  system elapsed 
# 8.105   0.125   4.351 
# openBlas
# user  system elapsed 
# 7.808   0.115   4.039 

# imac
# user  system elapsed 
# 19.444   0.085  19.617 
# accelerate
# user  system elapsed 
# 7.838   0.137   4.195 

system.time(res <- cd::spectrum_dispersion(cl, material, medium = medium))
# user  system elapsed 
# 18.081   0.247  18.442 
# vecLib
# user  system elapsed 
# 8.557   0.239   5.062 
# openBlas
# user  system elapsed 
# 8.705   0.250   5.064 

# imac
# user  system elapsed 
# 25.984   0.348  26.676 
# user  system elapsed 
# 14.446   0.397  10.855 

res1 <- cd:::cpp_dispersion_spectrum(kn, medium, cl$positions, Alpha, 
                                     cl$angles, 0.0, 0L, 
                                     0L, TRUE, FALSE, TRUE, 
                                     10, 1e-4, FALSE)

P <- cd:::cpp_polarization(E, AlphaBlocks)

cext <- cd:::cpp_extinction(kn, P, Ein)
cabs <- cd:::cpp_absorption(kn, P, E)

# Rprof(NULL)

subset(res, type=="cross-section")
res1[1:2]
cext[1]/N
cabs[1]/N

sessionInfo()
# imac
# R Under development (unstable) (2016-03-02 r70268)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.3 (El Capitan)
# 
# locale:
#   [1] en_NZ.UTF-8/en_NZ.UTF-8/en_NZ.UTF-8/C/en_NZ.UTF-8/en_NZ.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.3       reshape2_1.4.1   ggplot2_2.1.0    rgl_0.95.1441   
# [5] cd_0.0.1         dielectric_0.2.3
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.3      rngWELL_0.10-4   grid_3.3.0       gtable_0.2.0    
# [5] magrittr_1.5     scales_0.4.0     stringi_1.0-1    statmod_1.4.24  
# [9] randtoolbox_1.17 tools_3.3.0      stringr_1.0.0    munsell_0.4.3   
# [13] colorspace_1.2-6

