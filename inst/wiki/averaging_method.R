
## @knitr load
require(cda)
require(ggplot2)
require(microbenchmark)
require(xtable)

theme_set(theme_minimal())
display_benchmark <- function(x, unit = "t"){
    x$time <- microbenchmark:::convert_to_unit(x$time, unit)
    res <- aggregate(time ~ expr, x, function(.x) c(mean(.x), median(.x), min(.x), max(.x)))
    res <- cbind(res$expr, as.data.frame(res$time))
    colnames(res) <- c("expr",  "mean", "median", "min", "max")
    print(xtable(res), type = 'html', html.table.attributes = '')
}


## @knitr cluster
gold <- epsAu(seq(400, 900))

cl <- cluster_dimer(d=100, 
              dihedral=10*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)

# achiral cluster (plane of symmetry)
cl2 <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)



## @knitr comparison
params <- expand.grid(N=c(10, 100, 1000, 5000),
                       averaging=c("grid", "GL", "QMC"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl, material=gold)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(N), group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


## @knitr achiral
params <- expand.grid(N=c(10, 100, 1000, 5000),
                       averaging=c("grid", "GL", "QMC"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl2, material=gold)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(N), group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


