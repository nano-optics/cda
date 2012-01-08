library(cda)

wvl <- seq(400, 800)
gold <- epsAu(wvl)

onecluster <- function(N, Nq=30){
  clust <- makeSpheresCluster(N=N, radius=5e-3, R0=12e-3, pitch=15e-3,
                              delta=pi/2, right=TRUE)
  m <- circular_dichroism_spectrum(clust, gold, n=1.33, N=Nq)

  invisible(m)
}

params <- data.frame(N=seq(2,7))
comparison <- mdply(params, onecluster)

cs2mol <- (6.022e23 /2303 * 1e-8) # from micron^2 to cm^2

labs <- labs(y=expression(epsilon*" ("*10^3*M^-1*cm^-1*")"),
           x=expression(Energy*" (eV)"), linetype="handedness")

p <- 
ggplot(subset(comparison, variable=="extinction")) +
  facet_wrap(~type, scales="free") + labs + theme_minimal()+
  geom_path(aes(energy, value*cs2mol/1e3, colour=N, group=N))

p

ggsave("fig3.pdf",p)
