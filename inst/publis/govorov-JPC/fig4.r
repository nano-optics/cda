require(cda)

wvl <- seq(300, 800)
gold <- epsAu(wvl)

n <- sqrt(1.8)
rad <- 5e-3
step <- 25e-3
R0 <- 12e-3

clust0 <- makeSpheresCluster(N=4, radius=rad, R0=R0,
                             pitch=step, delta=pi/2, right=F)

Na <- 6.022e23 
fac <- (Na/2303 * 1e-8) # from micron^2 to cm^2

simulation <- function(N, ...){

  clust <- makeSpheresCluster(N=N, radius=rad, R0=R0,
                              pitch=step, delta=pi/2,
                              right=TRUE)

  m0 <- circular_dichroism_spectrum(clust, gold, n=n, N=30, ..., progress=FALSE)
  transform(m0, value=value*fac/N)
}

all <- mdply(data.frame(N=seq(4, 20)), simulation, .progress="text")
sall <- subset(all, type == "CD")

p <- 
ggplot(sall) + facet_grid(.~variable, scales="free")+
  geom_path(aes(wavelength, value,  colour=factor(N), linetype=p))+
  labs(y=expression(sigma[ext]*" (1/Mcm)"), colour="N",
           x=expression(energy*" /"*eV)) + theme_minimal() +
  scale_y_continuous(breaks=seq(-20000, 20000, by=5000))

p

ggsave("fig4.pdf", p)
