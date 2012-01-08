require(cda)

wvl <- seq(300, 800)
gold <- epsAu(wvl)

n <- sqrt(1.8)
rad <- 2.5e-3

Na <- 6.022e23 
fac <- (Na/2303 * 1e-8) # from micron^2 to cm^2

simulation <- function(remove=NULL){
clust <- makeSpheresCluster(N=10, radius=rad, R0=8.5e-3,
                            pitch=19.2e-3, delta=pi/3,
                            right=TRUE)
if(!is.null(remove)){
  clust$r <- clust$r[-remove,]
  clust$sizes <- clust$sizes[-remove,]
  clust$angles <- clust$angles[-remove,]
}

N <- NROW(clust$r)
m0 <- circular_dichroism_spectrum(clust, gold, n=n, N=30, progress=FALSE)
m0$remove <- paste("particles-", paste(remove, sep="-",collapse=""), sep="")
transform(m0, value=value*fac/N)
}

removal <- function(remove){
  all <- ldply(remove, simulation, .progress="text")
  subset(all, type == "CD" & variable == "extinction")
}
param.lists <- list(one =  list(NULL, 2, 3, 4, 5),
                    two =  list(NULL, c(2,3), c(2,5), c(2,7), c(3,6), c(3,8)),
                    three =  list(NULL, c(2,3,8), c(2,5,9), 4:6))

all.simulations <- ldply(param.lists, removal)
all.simulations$.id <- factor(all.simulations$.id, levels=c("one", "two", "three"))

p <- 
  ggplot(all.simulations) +
  facet_wrap(~.id, scales="free",ncol=1)+
  geom_path(aes(wavelength, value, colour=factor(remove)))+
  labs(y=expression(sigma[ext]*" (1/Mcm)"), colour="removed",
           x=expression(energy*" /"*eV)) + theme_minimal()
p

ggsave("fig8.pdf", p)
