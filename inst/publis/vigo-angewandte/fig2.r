library(cda)
library(ggplot2)
library(gridExtra)

gold <- epsAu(seq(400, 700))

govorov.params <- list(R0=12e-3, pitch=15e-3, N=5,
                             delta=pi/2, delta0=0, right=TRUE,
                             a=0.005, b=0.005, c=0.005, n=1.33, material=gold)


## to get the results of Govorov,  replace the scaling of CD by N to this (incorrect) factor

dubious <- 1.33^6

## define the cluster
cs2mol <- (6.022e23 /2303 * 1e-8) # from micron^2 to cm^2
onecluster <- function(N, ...){
  clust <- makeHelixCluster(N=N, ...)
  m <- circular_dichroism_spectrum(clust, gold, n=1.33, N=30)
  m$value <- m$value * cs2mol
  invisible(m)
}

params <- data.frame(N=seq(2,7))
comparison <- mdply(params, onecluster)
comparison2 <- mdply(params, onecluster, b=0.0045, c=0.0045)

abslabs <- labs(y=expression(epsilon*10^3*" (M)"),
           x=expression(Wavelength*" (nm)"), linetype="handedness")
cdlabs <-  labs(y=expression(CD/N*" ("*10^3*M^-1*cm^-1*")"), x=expression(Wavelength*" (nm)"), colour="N")

comparison$cluster <- "spheres"
comparison2$cluster <- "ellipsoids"

both <- rbind(comparison, comparison2)

ann.rods <- data.frame(cluster="ellipsoids", N=seq(2, 7),
                       lab= c("N=2", as.character(seq(3, 6)), "N=7"), 
                       x=c(560, 560, 580, 560, 580, 580),
                       y=c(0, 0.008,0.022, 0.014, 0.018,  0.027),
                       stringsAsFactors=FALSE)

ann.spheres <- data.frame(cluster="spheres", N=seq(2, 7),
                          lab= c("N=2, ", as.character(seq(3, 6)), "N=7"), 
                          x=c(550,563, 550, 550, 550, 550),
                          y=c(-0.02, -0.02, 0.12, -0.22, -0.1, 0.04)*1e-2,
                          stringsAsFactors=FALSE)

cs2mol <- (6.022e23 /2303 * 1e-8) # from micron^2 to cm^2

dummy <- data.frame(wavelength=550, CD=1.2*c(-6.2, 6, -62, 60), cluster=rep(c("spheres", "ellipsoids"), each=2))

p <- 
  ggplot(data=subset(both, variable=="extinction")) +
  facet_grid(type~cluster, scales="free_y")+
  scale_x_continuous(breaks=seq(0.4, 0.7, by=0.05)*1e3, expand=c(0, 0))+
  scale_y_continuous()+
  geom_line(aes(wavelength*1e3, value/N/1e3, colour=factor(N),group=N))+
  theme_minimal()+
  geom_blank(data=dummy, aes(wavelength, CD))+
  cdlabs + opts(strip.text.y=theme_blank(), strip.background=theme_blank(),
                panel.border = theme_rect(fill = NA, colour = "grey50"),
                panel.grid.major = theme_line(colour = "grey90", size = 0.2)) +
  opts(legend.position="none")

p

ggsave("fig2.pdf", p, width=7, height=7)


