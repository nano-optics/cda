library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 900))

original <- makeHelixCluster.zt(z=c(0, 0), theta=c(0, 10)*pi/180)
## copies
swapped <- shifted <- rotated <- original

##----------------------------------------------
## test invariance by swapping elements of the cluster
##----------------------------------------------

swapped$r <- original$r[2:1, ]
swapped$angles <- original$angles[2:1, ]

##----------------------------------------------
## test invariance by translation of the cluster
##----------------------------------------------

shifted$r <- original$r + 1

##--------------------------------------------
## test invariance by rotation of the cluster
##--------------------------------------------

rotated <- makeHelixCluster.zt(z=c(0, 0), theta=c(0, -10)*pi/180)

## model all clusters

tests <- llply(list(original=original, swapped=swapped, shifted=shifted, rotated=rotated),
             circular_dichroism_spectrum, material=gold, n=1.5, N=30,
                                   progress=FALSE, result.mat=FALSE)


m <- melt(tests, id=c("wavelength", "energy","variable","type"))
                                  
p <-
  ggplot(data=m)+ facet_grid(type~variable) +
  labs(y=expression(CD==bgroup("[", sigma[RH]-sigma[LH], "]")[Omega]*" /"*mu*m^2),
           x=expression(wavelength*" /"*mu*m))+ 
  geom_line(aes(wavelength, value, colour=L1))+ theme_bw()+
  labs(colour="cluster")

## ggsave("invariance.pdf", p)


all(sapply(tests, all.equal, target = tests[[1]]))
