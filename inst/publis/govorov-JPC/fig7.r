require(cda)

wvl <- seq(300, 800)
gold <- epsAu(wvl)


n <- sqrt(1.8)
rad <- 2.5e-3
## step <- 75e-3
R0 <- 8.5e-3

simulation <- function(step=7.5, ...){

  clust <- makeSpheresCluster(N=13, radius=rad, R0=R0,
                              pitch=step, delta=pi/3,
                              right=TRUE)

  Na <- 6.022e23 
  fac <- (Na/2303 * 1e-8) # from micron^2 to cm^2
  m0 <- circular_dichroism_spectrum(clust, gold, n=n, N=30, ..., progress=FALSE)
  transform(m0, value=value*fac/13)
}

st <- c(7.5,9,10.8,12.6,14.4,15.6,16.8,18,19.2,21,22.8,24.6,27,30,45,60)*1e-3
all <- mdply(data.frame(step=st), simulation, .progress="text")
sall <- subset(all, type == "CD" & variable == "extinction")

p <- 
ggplot(sall) + #facet_grid(.~variable, scales="free")+
  geom_path(aes(wavelength, value,  colour=factor(step*1e3), linetype=p))+
  labs(y=expression(sigma[ext]*" (1/Mcm)"), colour="pitch/nm",
           x=expression(energy*" /"*eV)) + theme_minimal()


p

ggsave("fig7a.pdf", p)

peak <- ddply(sall, .(step), summarize, peak = max(value))

p2 <- 
ggplot(peak,aes(step, peak)) +
  geom_path()+
  geom_point()+
  labs(x="pitch/nm", y="max(CD)") +
  theme_minimal()

ggsave("fig7b.pdf", p2)

